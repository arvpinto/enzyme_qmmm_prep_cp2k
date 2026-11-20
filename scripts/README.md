# scripts

| Script | Section | Description |
| :--- | :--- | :--- |
| prmtop_fix_types.py | [I](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_I.html) | renames the atom types whose CP2K kind names collide, or whose element CP2K would guess wrong (`-h` for options, `--dry-run` to only check) |
| mol2_qmsel.sh | [II](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_II.html) | QM selection for VMD, ParmEd and PyMOL from a mol2 file of the QM region |
| vmd_cp2k-qmmm.tcl | [III](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_III.html) | writes forceeval_subsys.cp2k.inc, forceeval_qmmm.cp2k.inc and qmmm-ee.parmed.in from a topology, coordinates and the QM selection |
| cp2k_funcs | [V](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_V.html) to [VIII](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_VIII.html) | bash functions to inspect CP2K outputs (`source cp2k_funcs; cp2k_help`) |
| vmd_cp2k-repr.tcl | [V](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_V.html) to [VIII](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_VIII.html) | loads a topology and CP2K DCD trajectories and draws the QM region with dynamic bonds and distance labels |
| plot_cp2k.py | [V](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_V.html) to [VII](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_VII.html) | energy against the iteration number of GEO_OPT runs (`opt`) and temperature against the time of MD runs (`temp`), shown and saved at 600 dpi (`-h` for options) |
| smd_jarzynski.py | [VIII](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_VIII.html) | work and Jarzynski free energy profile from steered MD metadynLog files (`-h` for options) |
