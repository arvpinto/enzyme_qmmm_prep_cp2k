# system_prep

System preparation of the tutorial (Sections I to IV): the force field topology, the coordinates, the QM region and the CP2K include files. Every CP2K run reads its topology and include files from this directory.

| File | Made by | Used for |
| :--- | :--- | :--- |
| frame_xray.pdb | starting structure (protein, oligosaccharide, Ca2+, Na+, Cl- and water) | `tleap -f leap_amyl.in` |
| leap_amyl.in | tleap input (ff14SB, GLYCAM_06j-1, TIP3P) | Section I |
| amylase_prep.prmtop, frame_xray.rst7 | `tleap -f leap_amyl.in` | topology and coordinates from tleap |
| amylase.prmtop | `python3 prmtop_fix_types.py -p amylase_prep.prmtop -o amylase.prmtop` | MM minimization and mechanical embedding QM/MM |
| type_renaming.dat | prmtop_fix_types.py | list of the atom types renamed for CP2K |
| hl_65atoms.mol2 | QM region (65 atoms), also shown in the 3D viewer of Section II | `./mol2_qmsel.sh -o qm_sel hl_65atoms.mol2` |
| qm_sel.vmd, qm_sel.mask, qm_sel.pml | mol2_qmsel.sh | QM selection for VMD, ParmEd (`netCharge`) and PyMOL |
| forceeval_subsys.cp2k.inc | vmd_cp2k-qmmm.tcl | &CELL of the system (MM minimization) |
| forceeval_qmmm.cp2k.inc | vmd_cp2k-qmmm.tcl | &CELL, &QM_KIND and &LINK sections of the QM/MM runs |
| qmmm-ee.parmed.in | vmd_cp2k-qmmm.tcl | charge correction for electrostatic embedding |
| qmmm-ee.parmed.log | `parmed -O -i qmmm-ee.parmed.in > qmmm-ee.parmed.log` | log of the correction |
| amylase_ee.prmtop | ParmEd | electrostatic embedding QM/MM runs |

Commands, in order (run in a working copy, with the scripts of [../scripts](../scripts)):

```
tleap -f leap_amyl.in
vmd -parm7 amylase_prep.prmtop -rst7 frame_xray.rst7
python3 prmtop_fix_types.py -p amylase_prep.prmtop -o amylase.prmtop
./mol2_qmsel.sh -o qm_sel hl_65atoms.mol2
echo "netCharge $(cat qm_sel.mask)" | parmed -n amylase.prmtop
vmd -parm7 amylase.prmtop -rst7 frame_xray.rst7 -e vmd_cp2k-qmmm.tcl -args -selfile qm_sel.vmd
parmed -O -i qmmm-ee.parmed.in > qmmm-ee.parmed.log
```

Tutorial: [I](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_I.html), [II](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_II.html), [III](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_III.html), [IV](https://arvpinto.github.io/enzyme_qmmm_prep_cp2k/enzyme_qmmm_prep_IV.html)
