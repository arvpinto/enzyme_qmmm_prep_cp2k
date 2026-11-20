# mm_em

MM minimization of the whole system with CP2K (METHOD FIST, 1000 L-BFGS steps). Tutorial: [Section V](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_V.html).

| File | Description |
| :--- | :--- |
| mm_em.inp | CP2K input |
| mm_em.out | CP2K output of the example (CP2K 2023.1, 64 MPI ranks) |
| mm_em_energy.dat, mm_em_energy.png | energy against the iteration number, from `python3 plot_cp2k.py opt mm_em.out` |
| MM_EM-1.restart | final structure, read by qmmm_em_me.inp (Section VI) |

Reads: amylase.prmtop, frame_xray.rst7 and forceeval_subsys.cp2k.inc (Sections I to III; also in ../system_prep).

In the tutorial every CP2K run is done in the same working directory (Section I), where all these files are, with a copy of the input:

```
cp ../enzyme_qmmm_prep_cp2k/mm_em/mm_em.inp .
cp2k.psmp -i mm_em.inp -o mm_em.out
```

For a parallel run, put the MPI launcher in front of cp2k.psmp (e.g. `mpirun -np 128 cp2k.psmp ...`).

The trajectory, MM_EM-pos-1.dcd, is attached to the [latest release](https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/releases/latest).
