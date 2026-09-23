# H2 + D2 -> 2 HD: isotope-aware thermochemistry

## Calculation

- Method: B3LYP/6-31G(d), PySCF 2.14; gas-phase RRHO treatment.
- Temperature / pressure: 298.15 K / 101325 Pa.
- Resources: one CPU core (`OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and `NUMEXPR_NUM_THREADS` set to 1).
- Electronic structure and Cartesian Hessian were optimized/evaluated once for the H--H Born--Oppenheimer surface. Nuclear masses were then set explicitly to 1.00782503223 u (H) and 2.01410177812 u (D) for isotope-specific translation, rotation, and vibration. Rotational symmetry numbers: H2=2, D2=2, HD=1.
- The optimized electronic energy is the same for all isotopologues: -1.1754824109 Eh. The optimized bond distance is 0.7428 Angstrom.

## Species values

Values are in kJ mol-1 for enthalpy and Gibbs components, and J mol-1 K-1 for entropy components. `elec` includes the Born--Oppenheimer electronic energy; `vib` includes ZPE and the finite-temperature vibrational term.

### Enthalpy components

| Species | elec | trans | rot | vib | total H |
| --- | ---: | ---: | ---: | ---: | ---: |
| H2 | -3086.228641 | 6.197391 | 2.478956 | 26.635632 | -3050.916662 |
| D2 | -3086.228641 | 6.197391 | 2.478956 | 18.841483 | -3058.710811 |
| HD | -3086.228641 | 6.197391 | 2.478956 | 23.070090 | -3054.482204 |

### Entropy components

| Species | elec | trans | rot | vib | total S |
| --- | ---: | ---: | ---: | ---: | ---: |
| H2 | 0.000000 | 117.487750 | 12.769431 | 0.000000 | 130.257181 |
| D2 | 0.000000 | 126.122883 | 18.526186 | 0.000034 | 144.649103 |
| HD | 0.000000 | 122.538197 | 20.922366 | 0.000001 | 143.460565 |

### Gibbs components

| Species | elec | trans | rot | vib | total G |
| --- | ---: | ---: | ---: | ---: | ---: |
| H2 | -3086.228641 | -28.831582 | -1.328250 | 26.635632 | -3089.752841 |
| D2 | -3086.228641 | -31.406147 | -3.044626 | 18.841473 | -3101.837941 |
| HD | -3086.228641 | -30.337373 | -3.759047 | 23.070089 | -3097.254972 |

## Reaction result: H2 + D2 -> 2 HD

| Quantity | elec | trans | rot | vib | total |
| --- | ---: | ---: | ---: | ---: | ---: |
| Delta H (kJ mol-1) | 0.000000 | 0.000000 | 0.000000 | 0.663064 | **0.663064** |
| Delta S (J mol-1 K-1) | 0.000000 | 1.465761 | 10.549115 | -0.000031 | **12.014845** |
| Delta G (kJ mol-1) | 0.000000 | -0.437017 | -3.145219 | 0.663073 | **-2.919162** |

Thus, at 298.15 K and 1 atm, the calculated reaction Gibbs free energy is **Delta G = -2.919 kJ mol-1**. The negative value is mainly driven by the rotational entropy gain of heteronuclear HD (rotational symmetry number 1 versus 2 for H2 and D2); the positive vibrational/ZPE enthalpy contribution partly offsets it.

## Comparison with the experimental equilibrium benchmark

For this gas-phase isotope exchange, the equilibrium constant is conventionally defined as `K = p(HD)^2 / [p(H2) p(D2)]`. A historical thermodynamic compilation used for hydrogen-isotope-mixture measurements reports **K = 3.26 at 298 K** for this reaction, compiled from Jones' 1949 data. [The compiled table](https://www.osti.gov/servlets/purl/4147038) identifies this reaction and value; a Sandia measurement-standard presentation independently lists the same 298 K constant and its `p(HD)^2/[p(H2)p(D2)]` definition. [Sandia HDT standards](https://www.osti.gov/servlets/purl/1420830)

Using `Delta G° = -RT ln K` at 298.15 K gives the following comparison:

| Quantity | This calculation | Literature benchmark | Difference (calculation - benchmark) |
| --- | ---: | ---: | ---: |
| K | 3.2465 | 3.2600 | -0.41% |
| Delta G° (kJ mol-1) | -2.9192 | -2.9295 | +0.0103 |

The calculated equilibrium free energy therefore agrees with the 298 K literature benchmark to within 0.011 kJ mol-1. This benchmark is a compiled thermodynamic equilibrium value rather than a new raw measurement in the cited source; its use here compares the predicted equilibrium composition directly with established isotope-exchange data.

## Files

- `run_isotope_thermochemistry.py`: reproducible one-core calculation script.
- `results.json`: unrounded machine-readable results, including harmonic frequencies (H2 4453.13, D2 3150.05, HD 3857.02 cm-1).

## MAESTRO support note

MAESTRO was used for capability triage. Its ThermoTask does not accept D as an isotope and does not expose isotope-aware thermochemistry; IsotopeShiftTask only produces isotope-shifted frequencies. The support gap is recorded in `../maestro-episodes/2026-09-23-h2-d2-hd-isotope-thermochemistry-gap.md`.
