# qmmm_md_eq

QM/MM MD of the reactant (electrostatic embedding, NVT, 300 K, 2000 steps of 1 fs). Tutorial: [Section VII](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_VII.html).

| File | Description |
| :--- | :--- |
| qmmm_md_eq.inp | CP2K input |
| qmmm_md_eq.out | CP2K output of the example (CP2K 2023.1, 128 MPI ranks) |
| QMMM_MD_EQ-1.ener | energies and temperature of every step |
| QMMM_MD_EQ-1.restart, QMMM_MD_EQ-1.restart.bin | final positions, velocities and thermostat, read by ../qmmm_md_smd/qmmm_md_smd.inp |
| QMMM_MD_EQ-RESTART.wfn | final wavefunction |

Needs: ../system_prep/amylase_ee.prmtop, ../system_prep/forceeval_qmmm.cp2k.inc, ../qmmm_em/QMMM_EM_EE-1.restart, ../qmmm_em/QMMM_EM_EE-RESTART.wfn

```
cp ../system_prep/{amylase_ee.prmtop,forceeval_qmmm.cp2k.inc} ../qmmm_em/{QMMM_EM_EE-1.restart,QMMM_EM_EE-RESTART.wfn} .
mpirun -np 128 cp2k.psmp -i qmmm_md_eq.inp -o qmmm_md_eq.out
```

The trajectory, QMMM_MD_EQ-pos-1.dcd, is attached to the [latest release](https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/releases/latest).
