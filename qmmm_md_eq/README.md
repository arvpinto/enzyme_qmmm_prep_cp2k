# qmmm_md_eq

QM/MM MD of the reactant (electrostatic embedding, NVT, 300 K, 2000 steps of 1 fs). Tutorial: [Section VII](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_VII.html).

| File | Description |
| :--- | :--- |
| qmmm_md_eq.inp | CP2K input |
| qmmm_md_eq.out | CP2K output of the example (CP2K 2023.1, 128 MPI ranks) |
| QMMM_MD_EQ-1.ener | energies and temperature of every step |
| QMMM_MD_EQ-1.restart, QMMM_MD_EQ-1.restart.bin | final positions, velocities and thermostat, read by qmmm_md_smd.inp (Section VIII) |
| QMMM_MD_EQ-RESTART.wfn | final wavefunction |
| qmmm_md_eq_temperature.dat, qmmm_md_eq_temperature.png | temperature of the system against time, from `python3 plot_cp2k.py temp QMMM_MD_EQ-1.ener --prefix qmmm_md_eq_temperature` |

Reads: amylase_ee.prmtop and forceeval_qmmm.cp2k.inc (Sections III and IV; also in ../system_prep), and QMMM_EM_EE-1.restart and QMMM_EM_EE-RESTART.wfn (Section VI).

In the tutorial every CP2K run is done in the same working directory (Section I), where all these files are, with a copy of the input:

```
cp ../enzyme_qmmm_prep_cp2k/qmmm_md_eq/qmmm_md_eq.inp .
cp2k.psmp -i qmmm_md_eq.inp -o qmmm_md_eq.out
```

For a parallel run, put the MPI launcher in front of cp2k.psmp (e.g. `mpirun -np 128 cp2k.psmp ...`).

The trajectory, QMMM_MD_EQ-pos-1.dcd, is attached to the [latest release](https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/releases/latest).
