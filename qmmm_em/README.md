# qmmm_em

QM/MM minimizations (PBE/DZVP-GTH-PBE with GPW, 65 QM atoms). Tutorial: [Section VI](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_VI.html).

| File | Description |
| :--- | :--- |
| qmmm_em_me.inp, qmmm_em_me.out | mechanical embedding (ECOUPL NONE), 500 steps, from MM_EM-1.restart |
| QMMM_EM_ME-1.restart, QMMM_EM_ME-RESTART.wfn | final structure and wavefunction, read by qmmm_em_ee.inp |
| qmmm_em_ee.inp, qmmm_em_ee.out | electrostatic embedding (ECOUPL GAUSS), converged |
| qmmm_em_energy.dat, qmmm_em_energy.png | energy against the iteration number of both runs, from `python3 plot_cp2k.py opt qmmm_em_me.out qmmm_em_ee.out --prefix qmmm_em_energy` |
| QMMM_EM_EE-1.restart, QMMM_EM_EE-RESTART.wfn | final structure and wavefunction, read by qmmm_md_eq.inp (Section VII) |

Reads: amylase.prmtop (ME), amylase_ee.prmtop (EE) and forceeval_qmmm.cp2k.inc (Sections I to IV; also in ../system_prep), and MM_EM-1.restart (Section V).

In the tutorial every CP2K run is done in the same working directory (Section I), where all these files are, with a copy of the input:

```
cp ../enzyme_qmmm_prep_cp2k/qmmm_em/{qmmm_em_me.inp,qmmm_em_ee.inp} .
cp2k.psmp -i qmmm_em_me.inp -o qmmm_em_me.out
cp2k.psmp -i qmmm_em_ee.inp -o qmmm_em_ee.out
```

For a parallel run, put the MPI launcher in front of cp2k.psmp (e.g. `mpirun -np 128 cp2k.psmp ...`).

The trajectories, QMMM_EM_ME-pos-1.dcd and QMMM_EM_EE-pos-1.dcd, are attached to the [latest release](https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/releases/latest).
