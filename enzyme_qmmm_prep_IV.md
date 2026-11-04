<br>
<h2> <p align="center"> <b>IV - Charge correction for electrostatic embedding</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
amylase.prmtop  ->  qmmm-ee.parmed.in (ParmEd)  ->  amylase_ee.prmtop  ->  &QMMM ECOUPL GAUSS
</pre>

---

<br/>

<h3> <p align="center"><b> Why the MM charges are corrected </b></p> </h3>

<br/>

- <p align="justify">With mechanical embedding (ECOUPL NONE, used in the first QM/MM minimization), the QM and MM regions interact through the force field: the QM atoms keep their MM point charges for the QM-MM terms. With electrostatic embedding (ECOUPL GAUSS, used in every other QM/MM step), the QM electron density feels the MM point charges directly and the MM charges of the QM atoms are removed: the QM atoms are represented by their nuclei and their electrons.</p>

- <p align="justify">Residues have integer charges in the force field, but their fragments do not. When the boundary cuts a residue, the part left in the MM region carries a fractional charge, right next to the QM density. The correction restores an integer charge in the MM part of every residue cut by the boundary, by spreading the residual evenly over its MM atoms, as in the <a href="https://manual.cp2k.org/trunk/methods/qm_mm/builtin.html" target="_blank">CP2K biochemical QM/MM example</a>.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> The correction of the example </b></p> </h3>

<br/>

<p align="justify">The last part of the output of vmd_cp2k-qmmm.tcl (Section III) reports the MM charge of the selection, -2.2063 e, and the correction of each residue that is cut:</p>

| Residue | MM charge of the QM part | Integer | MM atoms of the residue | Correction per MM atom |
| :--- | ---: | ---: | :--- | ---: |
| ASP196 | -0.8581 | -1 | N, H, CA, HA, C, O (6) | +0.02365 |
| GLH232 | -0.0221 | 0 | N, H, CA, HA, CB, HB2, HB3, C, O (9) | -0.00246 |
| ASP299 | -0.8581 | -1 | N, H, CA, HA, C, O (6) | +0.02365 |
| 4GA499 | -0.4680 | 0 | all atoms except O4 (20) | -0.02340 |

<br/>

<p align="justify">4GA500, 4GA501 and WAT511 are whole in the QM region, so they need no correction. After the correction the MM part of each cut residue has the integer charge of the residue minus that of its QM part (0 for all four), the classical system has a charge of +2 and the QM region carries -2, the &DFT CHARGE of the inputs.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> The corrected topology </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
cat qmmm-ee.parmed.in
parm amylase.prmtop

change charge @3018 -0.4926498184601466
change charge @3019 0.31724992891152698
change charge @3020 0.061750018348296479
...
change charge @7677 -0.59739987403154371
changeLJSingleType @%HO,ho 0.3019 0.047
changeLJSingleType @%HW,hw 0.3019 0.047

outparm amylase_ee.prmtop
</pre>

<pre style="color: white; background-color: black;">
parmed -O -i qmmm-ee.parmed.in > qmmm-ee.parmed.log
</pre>

<pre style="color: white; background-color: black;">
cat qmmm-ee.parmed.log
...
Adding prmtop amylase.prmtop to parm list. amylase.prmtop is the active parm.
Changing CHARGE of atom # 3018 (N) from -0.5163 to -0.4926
Changing CHARGE of atom # 3019 (H) from 0.2936 to 0.3172
...
Changing CHARGE of atom # 7677 (O5) from -0.5740000000000001 to -0.5974
Changing @%HO,ho Lennard-Jones well depth from 0.0000 to 0.0470 (kal/mol) and radius from 0.0000 to 0.3019 (Angstroms)
Changing @%HW,hw Lennard-Jones well depth from 0.0000 to 0.0470 (kal/mol) and radius from 0.0000 to 0.3019 (Angstroms)
Outputting Amber topology file amylase_ee.prmtop
Done!
</pre>

<br/>

- <p align="justify">The script has one change charge line per MM atom of the cut residues (41 lines here, by serial number) and two changeLJSingleType lines. The charges are written with full precision; the small differences from the four-decimal values of the force field come from the single-precision charges that VMD stores.</p>

- <p align="justify">The changeLJSingleType lines give a Lennard-Jones radius and well depth (R<sub>min</sub>/2 = 0.3019 &Aring;, &epsilon; = 0.047 kcal/mol, the values of the CP2K example) to the hydroxyl hydrogens of the protein (type HO) and to the TIP3P hydrogens (HW), which have none in the force field. Without repulsion, these positive MM charges can get too close to the QM density in electrostatic embedding QM/MM MD. The change applies to every atom of those types, including all water molecules, so it modifies the MM model: use it consistently in all the calculations that you compare. The GLYCAM hydroxyl hydrogens (Ho, renamed HT in Section I) already have LJ parameters in GLYCAM. Change the list with -lj in vmd_cp2k-qmmm.tcl, and note that the masks are case sensitive and use the type names of the topology the script reads.</p>

- <p align="justify">amylase_ee.prmtop is used only for electrostatic embedding, as both PARM_FILE_NAME and CONN_FILE_NAME (the inputs set both from @SET PARM). The MM minimization and the mechanical embedding minimization keep amylase.prmtop. The QM atoms keep their original charges in amylase_ee.prmtop, so its total charge is not an integer (-0.2063); that is expected, since CP2K removes them.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Checks </b></p> </h3>

<br/>

- <p align="justify">The correction creates new charges, and with them new CP2K kinds. Check the corrected topology with prmtop_fix_types.py in dry-run mode: with the types renamed in Section I, there is nothing left to fix.</p>

<pre style="color: white; background-color: black;">
python3 prmtop_fix_types.py -p amylase_ee.prmtop --dry-run -q
</pre>

<pre style="color: white; background-color: black;">
prmtop_fix_types.py v2.0
amylase_ee.prmtop: 43 atom types, 284 CP2K kinds
  kind names shared by different types: none

Types whose element CP2K cannot guess (give them a &amp;KIND section with ELEMENT in &amp;SUBSYS):
  Cl-    element Cl  kind name(s): CL-
  Na+    element Na  kind name(s): NA+

No atom type needs to be renamed.
</pre>

<br/>

- <p align="justify">In the output of the electrostatic embedding runs (Sections VI to VIII), CP2K prints the charge of the classical system once the QM atoms are removed. It has to be an integer, +2 here; with the uncorrected topology it would be +2.2063.</p>

<pre style="color: white; background-color: black;">
grep "Total Charge of the Classical System" qmmm_em_ee.out
 CHARGE_INFO| Total Charge of the Classical System:                     2.000000
</pre>

<br/>
