---
task: unsupported
engine: none
error_class: support_gap
outcome: support_gap
---
Symptom: The requested reaction H2 + D2 -> 2 HD needs isotope-aware Gibbs thermochemistry. ThermoTask rejected XYZ element label D during charge/spin suggestion.

Attempts: Queried ThermoTask, IsotopeShiftTask, and ThermoReevalTask. ThermoTask produces Gibbs free energy and entropy components but accepts only ordinary element labels. IsotopeShiftTask supports mass-number substitutions only for harmonic frequencies. ThermoReevalTask accepts frequencies, geometry, and electronic energy but exposes no isotope-substitution or mass input.

Result: No registered MAESTRO task combines isotopic masses with translational, rotational, vibrational, and electronic thermochemistry.

Context: Direct PySCF workflow is required for the requested H2/D2/HD reaction free energy at the selected DFT B3LYP/6-31G(d) level.
