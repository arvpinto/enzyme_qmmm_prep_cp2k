# qmmm_md_smd

Steered QM/MM MD along CV = d(O4-C1) - d(HE2(Glu232)-O4), from -2 to +2 angstrom in 2 ps, and its Jarzynski work profile. Tutorial: [Section VIII](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_VIII.html).

| File | Description |
| :--- | :--- |
| qmmm_md_smd.inp | CP2K input |
| qmmm_md_smd.out | CP2K output of the example (CP2K 2023.1, 128 MPI ranks) |
| QMMM_MD_SMD-COLVAR.metadynLog | time (fs) and CV (bohr) every 10 steps |
| QMMM_MD_SMD-1.ener | energies and temperature of every step |
| QMMM_MD_SMD-1.restart, QMMM_MD_SMD-1.restart.bin, QMMM_MD_SMD-RESTART.wfn | final state (glycosyl-enzyme intermediate) |
| smd_jarzynski.dat, smd_jarzynski.png | work profile from smd_jarzynski.py |

Reads: amylase_ee.prmtop and forceeval_qmmm.cp2k.inc (Sections III and IV; also in ../system_prep), and QMMM_MD_EQ-1.restart, QMMM_MD_EQ-1.restart.bin and QMMM_MD_EQ-RESTART.wfn (Section VII).

In the tutorial every CP2K run is done in the same working directory (Section I), where all these files are, with a copy of the input:

```
cp ../enzyme_qmmm_prep_cp2k/qmmm_md_smd/qmmm_md_smd.inp .
cp2k.psmp -i qmmm_md_smd.inp -o qmmm_md_smd.out
python3 smd_jarzynski.py -l QMMM_MD_SMD-COLVAR.metadynLog -o qmmm_md_smd.out
```

For a parallel run, put the MPI launcher in front of cp2k.psmp (e.g. `mpirun -np 128 cp2k.psmp ...`).

The trajectory, QMMM_MD_SMD-pos-1.dcd, is attached to the [latest release](https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/releases/latest).
