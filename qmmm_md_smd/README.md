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

Needs: ../system_prep/amylase_ee.prmtop, ../system_prep/forceeval_qmmm.cp2k.inc, ../qmmm_md_eq/QMMM_MD_EQ-1.restart, ../qmmm_md_eq/QMMM_MD_EQ-1.restart.bin, ../qmmm_md_eq/QMMM_MD_EQ-RESTART.wfn

```
cp ../system_prep/{amylase_ee.prmtop,forceeval_qmmm.cp2k.inc} .
cp ../qmmm_md_eq/{QMMM_MD_EQ-1.restart,QMMM_MD_EQ-1.restart.bin,QMMM_MD_EQ-RESTART.wfn} .
mpirun -np 128 cp2k.psmp -i qmmm_md_smd.inp -o qmmm_md_smd.out
python3 ../scripts/smd_jarzynski.py -l QMMM_MD_SMD-COLVAR.metadynLog -o qmmm_md_smd.out
```

The trajectory, QMMM_MD_SMD-pos-1.dcd, is attached to the [latest release](https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/releases/latest).
