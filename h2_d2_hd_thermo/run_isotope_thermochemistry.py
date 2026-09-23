"""Isotope-aware B3LYP/6-31G(d) thermochemistry for H2 + D2 -> 2 HD."""

import json
import os
from pathlib import Path

os.environ.update({
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
})

import numpy as np
from pyscf import dft, gto
from pyscf.data import nist
from pyscf.geomopt.geometric_solver import optimize
from pyscf.hessian import thermo


ROOT = Path(__file__).resolve().parent
TEMPERATURE = 298.15
PRESSURE_PA = 101325.0
H_MASS = 1.00782503223
D_MASS = 2.01410177812
HARTREE_TO_KJMOL = nist.HARTREE2J * nist.AVOGADRO / 1000.0
EH_PER_K_TO_JMOLK = nist.HARTREE2J * nist.AVOGADRO


def make_mf(mol):
    return dft.RKS(mol).set(xc="B3LYP", conv_tol=1e-10)


def isotope_thermo(mol, energy, hessian, masses, symmetry_number):
    """RRHO components with specified nuclear masses and rotational symmetry."""
    freq = thermo.harmonic_analysis(mol, hessian, mass=np.array(masses))["freq_au"]
    coords = mol.atom_coords()
    mass = np.array(masses)
    center = np.einsum("z,zx->x", mass, coords) / mass.sum()
    coords = coords - center

    k_b, h = nist.BOLTZMANN, nist.PLANCK
    r_eh = k_b / nist.HARTREE2J
    components_h = {"elec": energy}
    components_s = {"elec": r_eh * np.log(mol.multiplicity)}

    total_mass = mass.sum() * nist.ATOMIC_MASS
    q_trans = ((2 * np.pi * total_mass * k_b * TEMPERATURE / h**2) ** 1.5
               * k_b * TEMPERATURE / PRESSURE_PA)
    components_h["trans"] = 2.5 * r_eh * TEMPERATURE
    components_s["trans"] = r_eh * (2.5 + np.log(q_trans))

    rot_const = thermo.rotation_const(mass, coords, "GHz")
    b_hz = rot_const[1] * 1e9
    q_rot = k_b * TEMPERATURE / (symmetry_number * h * b_hz)
    components_h["rot"] = r_eh * TEMPERATURE
    components_s["rot"] = r_eh * (1 + np.log(q_rot))

    au2hz = (nist.HARTREE2J / (nist.ATOMIC_MASS * nist.BOHR_SI**2)) ** 0.5 / (2 * np.pi)
    vib_temp = freq.real[freq.real > 0] * au2hz * h / k_b
    reduced_temp = vib_temp / TEMPERATURE
    boltzmann = np.exp(-reduced_temp)
    zpe = 0.5 * r_eh * vib_temp.sum()
    components_h["vib"] = zpe + r_eh * TEMPERATURE * (reduced_temp * boltzmann / (1 - boltzmann)).sum()
    components_s["vib"] = r_eh * (reduced_temp * boltzmann / (1 - boltzmann) - np.log(1 - boltzmann)).sum()

    components_g = {key: components_h[key] - TEMPERATURE * components_s[key]
                    for key in components_h}
    return {
        "frequency_cm-1": (freq.real * au2hz / nist.LIGHT_SPEED_SI * 1e-2).tolist(),
        "zpe_hartree": float(zpe),
        "H_components_hartree": {key: float(value) for key, value in components_h.items()},
        "S_components_Eh_per_K": {key: float(value) for key, value in components_s.items()},
        "G_components_hartree": {key: float(value) for key, value in components_g.items()},
        "H_total_hartree": float(sum(components_h.values())),
        "S_total_Eh_per_K": float(sum(components_s.values())),
        "G_total_hartree": float(sum(components_g.values())),
    }


def reaction_value(results, key):
    return 2 * results["HD"][key] - results["H2"][key] - results["D2"][key]


def main():
    mol = gto.M(atom="H 0 0 0; H 0 0 0.740", basis="6-31G(d)", charge=0, spin=0)
    mol_opt = optimize(make_mf(mol), maxsteps=50)
    mf = make_mf(mol_opt)
    energy = mf.kernel()
    hessian = mf.Hessian().kernel()

    specs = {"H2": ([H_MASS, H_MASS], 2), "D2": ([D_MASS, D_MASS], 2), "HD": ([H_MASS, D_MASS], 1)}
    results = {name: isotope_thermo(mol_opt, energy, hessian, masses, sigma)
               for name, (masses, sigma) in specs.items()}
    delta = {
        "H_hartree": reaction_value(results, "H_total_hartree"),
        "S_Eh_per_K": reaction_value(results, "S_total_Eh_per_K"),
        "G_hartree": reaction_value(results, "G_total_hartree"),
        "H_kJ_mol": reaction_value(results, "H_total_hartree") * HARTREE_TO_KJMOL,
        "S_J_mol_K": reaction_value(results, "S_total_Eh_per_K") * EH_PER_K_TO_JMOLK,
        "G_kJ_mol": reaction_value(results, "G_total_hartree") * HARTREE_TO_KJMOL,
    }
    payload = {"method": "B3LYP/6-31G(d)", "temperature_K": TEMPERATURE,
               "pressure_Pa": PRESSURE_PA, "optimized_geometry_bohr": mol_opt.atom_coords().tolist(),
               "species": results, "reaction": delta}
    (ROOT / "results.json").write_text(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
