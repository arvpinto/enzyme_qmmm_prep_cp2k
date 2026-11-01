<br>
<h2> <p align="center"> <b>III - Include files for the &QMMM and &SUBSYS sections</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
&FORCE_EVAL
  METHOD QMMM
  ...
  &QMMM
    ...
    @INCLUDE forceeval_qmmm.cp2k.inc     ! QM cell, &QM_KIND, &LINK
  &END QMMM
  &SUBSYS
    @INCLUDE forceeval_subsys.cp2k.inc   ! MM cell (MM minimization)
...
</pre>

---

<br/>

<p align="justify">The &QMMM section tells CP2K which atoms are QM (&QM_KIND), how the bonds cut by the boundary are capped (&LINK) and the size of the box in which the QM problem is solved (&CELL). With 65 QM atoms these are long lists of atom numbers, so they are written by a script and included in every input with @INCLUDE. vmd_cp2k-qmmm.tcl writes them from the VMD selection of Section II, together with the MM cell and the ParmEd script of the charge correction (Section IV).</p>

<br/>

<h3> <p align="center"><b> Running the script </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
vmd -parm7 amylase.prmtop -rst7 frame_xray.inpcrd -e vmd_cp2k-qmmm.tcl -args -selfile qm_sel.vmd
</pre>

<p align="justify">or without arguments, typing (or pasting) the selection in the terminal when the script asks for it:</p>

<pre style="color: white; background-color: black;">
vmd -parm7 amylase.prmtop -rst7 frame_xray.inpcrd -e vmd_cp2k-qmmm.tcl

1) MAKE SURE YOUR SYSTEM IS THE TOP MOLECULE!
CONTROL STRING: TOP MOLECULE IS 0!

2) INTRODUCE SELECTION FOR QM LAYER IN VMD SELECTION FORMAT:
(name CB CG OD1 OD2 HB2 HB3 and resname ASP and resid 196) or ... or (name O H1 H2 and resname WAT and resid 511)
CONTROL STRING: 65 atoms were selected!
CONTROL STRING: 5 QM/MM bonds were cut (one &amp;LINK section each)!

3) DOUBLE-CHECK THAT INTEGER CHARGE OF SELECTION MATCHES ITS PARCELS!
...
</pre>

<br/>

- <p align="justify">The options (after -args) are -sel or -selfile for the selection, -pad for the padding of the QM cell on each side (default 5 &Aring;), -scale for QMMM_SCALE_FACTOR (default 0.5), -lj for the hydrogen types that get Lennard-Jones parameters in the ParmEd script (Section IV), -overwrite to replace existing output files (by default the script stops rather than overwrite them) and -quit to close VMD at the end. Add -dispdev text to run it without graphics.</p>

- <p align="justify">The coordinates must have the box information (tleap writes it at the end of the inpcrd file), since the MM cell is taken from them. In graphics mode, the script also draws the QM layer, the residues it cuts and the QM/MM bonds, with a label on each cut bond: rotate the active site and check the selection before going further.</p>

<br/>

<p align="justify">Three files are written (the versions of the example are in the example folder of <a href="CP2K_INPUTS.zip">CP2K_INPUTS.zip</a>):</p>

| File | Content | Used in |
| :--- | :--- | :--- |
| forceeval_subsys.cp2k.inc | &CELL of the MM system (the box of the coordinates) | em.inp (&SUBSYS) |
| forceeval_qmmm.cp2k.inc | QM &CELL, one &QM_KIND per element and one &LINK per cut bond | every QM/MM input (&QMMM) |
| qmmm-ee.parmed.in | charge correction for electrostatic embedding | Section IV |

<br/>

---

<br/>

<h3> <p align="center"><b> QM atoms </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
! CP2K > FORCE_EVAL > QMMM > QM_KIND
&QM_KIND H   !QM_KIND guess from element mass 1.0080000162124634
  MM_INDEX 3023 3024 3602 3603 3607 4651 4652 7680 7682 7684 7686 7688 7690 7692 7694 7695 7697 7701 7703 7705 7707 7709 7711 7713 7715 7716 7718 7784 7785
&END QM_KIND
&QM_KIND C   !QM_KIND guess from element mass 12.010000228881836
  MM_INDEX 3022 3025 3601 3604 4650 4653 7679 7681 7685 7689 7691 7693 7700 7702 7706 7710 7712 7714
&END QM_KIND
&QM_KIND O   !QM_KIND guess from element mass 16.0
  MM_INDEX 3026 3027 3605 3606 4654 4655 7678 7683 7687 7696 7698 7699 7704 7708 7717 7719 7720 7783
&END QM_KIND
</pre>

<br/>

- <p align="justify">MM_INDEX lists the serial numbers (1-based) of the QM atoms in the topology. The element is guessed from the mass, which is safer than the atom type: the QM kinds are named C, H and O, and these are the names of the &KIND sections that give the basis sets and pseudopotentials in the QM/MM inputs (DZVP-GTH-PBE for C and O, DZV-GTH-PBE for H). A QM region with other elements needs one &KIND section per element.</p>

- <p align="justify">The serials of the QM atoms can be read back from the include file with the cp2k_qmidx function of cp2k_funcs, which gives a VMD selection (used in Section VIII):</p>

<pre style="color: white; background-color: black;">
source cp2k_funcs
cp2k_qmidx forceeval_qmmm.cp2k.inc
serial 3023 3024 3602 3603 3607 4651 4652 7680 ... 7719 7720 7783
</pre>

<br/>

---

<br/>

<h3> <p align="center"><b> QM/MM links </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
! CP2K > FORCE_EVAL > QMMM > LINK
&LINK
  ALPHA_IMOMM 1.40     !CT-CT bond (default is ALPHA_IMOMM 1.38 for C-C bond)
  LINK_TYPE IMOMM
  MM_INDEX  3020
  QM_INDEX  3022
  QMMM_SCALE_FACTOR 0.5
&END LINK
...
&LINK
  ALPHA_IMOMM 1.44     !Cg-Os bond (default is ALPHA_IMOMM 1.38 for C-C bond)
  LINK_TYPE IMOMM
  MM_INDEX  7668
  QM_INDEX  7678
  QMMM_SCALE_FACTOR 0.5
&END LINK
&LINK
  ALPHA_IMOMM 1.45     !Os-Cg bond (default is ALPHA_IMOMM 1.38 for C-C bond)
  LINK_TYPE IMOMM
  MM_INDEX  7721
  QM_INDEX  7720
  QMMM_SCALE_FACTOR 0.5
&END LINK
</pre>

<br/>

| QM atom | MM atom | Bond | ALPHA_IMOMM |
| :--- | :--- | :--- | ---: |
| ASP196 CB (3022) | ASP196 CA (3020) | CT-CT | 1.40 |
| GLH232 CG (3601) | GLH232 CB (3598) | CT-CT | 1.39 |
| ASP299 CB (4650) | ASP299 CA (4648) | CT-CT | 1.40 |
| 4GA499 O4 (7678) | 4GA499 C4 (7668) | Os-Cg | 1.44 |
| 4GA501 O4 (7720) | 4GA502 C1 (7721) | Os-Cg | 1.45 |

<br/>

- <p align="justify">Each QM atom bonded to an MM atom is capped with a hydrogen placed along the QM-MM bond, at the QM-MM distance divided by ALPHA_IMOMM, and the forces on the capping hydrogen are projected onto the QM and MM atoms (IMOMM scheme). The CP2K default (1.38) corresponds to a C-C bond capped by a C-H bond. The script estimates ALPHA_IMOMM from the length of each cut bond in the input coordinates, as r / (1.090 - (1.526 - r)), a rough guess built for C-C cuts. For the C-C cuts it gives 1.39-1.40 (1.526/1.090 = 1.40 from the force field distances). For the two C-O cuts it gives 1.44-1.45, which puts the capping hydrogen of the boundary oxygen at about 0.98 &Aring;; the ratio of the GLYCAM equilibrium distances, r<sub>0</sub>(Cg-Os)/r<sub>0</sub>(Oh-Ho) = 1.46/0.96 = 1.52, is the value that reproduces the O-H equilibrium distance when the C-O bond is at its own equilibrium length. Check the capping O-H distances in the QM/MM minimization (Section VI) and adjust ALPHA_IMOMM if needed.</p>

- <p align="justify">QMMM_SCALE_FACTOR scales the charge of the MM link atom in the QM/MM electrostatic potential only. The MM link atom is about 1.5 &Aring; from a QM atom and its full point charge can overpolarize the QM density, so it is halved here (-scale changes it for all links). Keep the value fixed for all the calculations that you will compare.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> QM and MM cells </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
! CP2K > FORCE_EVAL > QMMM > CELL
&CELL
  ABC [angstrom] 19.038997650146484  19.068000793457031  20.088996887207031     !size of cell
  ALPHA_BETA_GAMMA 90.0 90.0 90.0      !shape of cell
&END CELL
</pre>

<pre style="color: white; background-color: black;">
! CP2K > FORCE_EVAL > SUBSYS > CELL
&CELL
  ABC [angstrom] 88.588600  102.292603  90.199600     !size of cell
  ALPHA_BETA_GAMMA 90.000000 90.000000 90.000000      !shape of cell
&END CELL
</pre>

<br/>

- <p align="justify">The QM cell is the extent of the QM atoms (about 9 x 9 x 10 &Aring;) plus 5 &Aring; on each side. CP2K centers the QM atoms in it at every step (CENTER EVERY_STEP in the inputs) and warns if a QM atom gets too close to its border. With GPW the electron density is expanded in plane waves in this cell, so the cell sets most of the cost of the QM part (the number of grid points grows with its volume) and must hold the density with some margin. Test the energy against the padding (-pad) before production runs, and remember that the sugars and the water can move during MD: if the QM region stretches, the cell has to grow with it.</p>

- <p align="justify">The QM region is charged (-2), and the inputs solve the QM electrostatics with the default periodic Poisson solver in the QM cell, without the &QMMM/&PERIODIC section (the QM/MM electrostatics are computed in real space with the GEEP Gaussian expansion of the MM charges). Images of a charged QM region interact through the periodic QM cell, so a small padding can bias the energies; the size test above also checks this. &QMMM/&PERIODIC (decoupling and recoupling of the QM images) is the alternative described in the CP2K manual.</p>

- <p align="justify">The MM cell is the box of the coordinates (setBox vdw in Section I). It is only needed in em.inp: the QM/MM inputs read the cell and the positions from the restart file of the previous step.</p>

<br/>
