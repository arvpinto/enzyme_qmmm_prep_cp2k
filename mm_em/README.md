# mm_em

MM minimization of the whole system with CP2K (METHOD FIST, 1000 L-BFGS steps). Tutorial: [Section V](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_V.html).

| File | Description |
| :--- | :--- |
| mm_em.inp | CP2K input |
| mm_em.out | CP2K output of the example (CP2K 2023.1, 64 MPI ranks) |
| MM_EM-1.restart | final structure, read by ../qmmm_em/qmmm_em_me.inp |

Needs, from ../system_prep: amylase.prmtop, frame_xray.rst7, forceeval_subsys.cp2k.inc

```
cp ../system_prep/{amylase.prmtop,frame_xray.rst7,forceeval_subsys.cp2k.inc} .
mpirun -np 64 cp2k.psmp -i mm_em.inp -o mm_em.out
```

The trajectory, MM_EM-pos-1.dcd, is attached to the [latest release](https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/releases/latest).
