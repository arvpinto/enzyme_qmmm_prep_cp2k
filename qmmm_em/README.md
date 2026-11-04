# qmmm_em

QM/MM minimizations (PBE/DZVP-GTH-PBE with GPW, 65 QM atoms). Tutorial: [Section VI](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_VI.html).

| File | Description |
| :--- | :--- |
| qmmm_em_me.inp, qmmm_em_me.out | mechanical embedding (ECOUPL NONE), 500 steps, from ../mm_em/MM_EM-1.restart |
| QMMM_EM_ME-1.restart, QMMM_EM_ME-RESTART.wfn | final structure and wavefunction, read by qmmm_em_ee.inp |
| qmmm_em_ee.inp, qmmm_em_ee.out | electrostatic embedding (ECOUPL GAUSS), converged |
| QMMM_EM_EE-1.restart, QMMM_EM_EE-RESTART.wfn | final structure and wavefunction, read by ../qmmm_md_eq/qmmm_md_eq.inp |

Needs, from ../system_prep: amylase.prmtop (ME), amylase_ee.prmtop (EE), forceeval_qmmm.cp2k.inc; and ../mm_em/MM_EM-1.restart

```
cp ../system_prep/{amylase.prmtop,amylase_ee.prmtop,forceeval_qmmm.cp2k.inc} ../mm_em/MM_EM-1.restart .
mpirun -np 128 cp2k.psmp -i qmmm_em_me.inp -o qmmm_em_me.out
mpirun -np 128 cp2k.psmp -i qmmm_em_ee.inp -o qmmm_em_ee.out
```

The trajectories, QMMM_EM_ME-pos-1.dcd and QMMM_EM_EE-pos-1.dcd, are attached to the [latest release](https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/releases/latest).
