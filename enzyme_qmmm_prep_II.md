<br>
<h2> <p align="center"> <b>II - Selection of the QM region</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
amylase.prmtop + frame_xray.inpcrd  ->  QM region picked in PyMOL (hl_65atoms.mol2)  ->  mol2_qmsel.sh  ->  VMD / PyMOL / Amber selections
</pre>

---

<br/>

<p align="justify">The QM region decides which chemistry the model can describe, how the environment polarizes the reacting atoms and how much each MD step costs. There is no universal recipe, but there are well established rules.</p>

<br/>

<h3> <p align="center"><b> Guidelines </b></p> </h3>

<br/>

- <p align="justify">Include every atom that forms or breaks a bond or changes its protonation state along the reaction, together with the groups conjugated with them (carboxylates, amides, aromatic rings), and the groups that interact strongly with the reacting atoms and change along the reaction (hydrogen bond donors to a developing charge, metal ions and their ligands, catalytic waters).</p>

- <p align="justify">Cut covalent bonds preferably between apolar sp<sup>3</sup> carbon atoms (C&alpha;-C&beta; or C&beta;-C&gamma; of side chains), never in a peptide bond, and keep the boundary at least three bonds away from any bond that changes during the reaction. Avoid a QM atom with two MM neighbours and an MM atom bonded to two QM atoms.</p>

- <p align="justify">Water molecules can be part of the QM region, but in QM/MM MD they can diffuse away. Check where they are at the end of each run.</p>

- <p align="justify">The cost of QM/MM MD grows quickly with the number of QM atoms (and with the size of the QM cell, Section III). A smaller, well chosen region with more sampling is often better than a larger one with less.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> QM region of the example </b></p> </h3>

<br/>

<p align="justify">The QM region of the example has 65 atoms (18 C, 18 O and 29 H) and a charge of -2:</p>

| Residue | QM atoms | Role |
| :--- | :--- | :--- |
| ASP196 | side chain from C&beta; (6 atoms) | nucleophile (Asp197 in HPA), 3.1 &Aring; from C1 of the -1 sugar |
| GLH232 | side chain from C&gamma; (7 atoms) | general acid: its HE2 is 2.0 &Aring; from the glycosidic O4 |
| ASP299 | side chain from C&beta; (6 atoms) | binds the O2 and O3 hydroxyl groups of the -1 sugar |
| 4GA499 | O4 only | glycosidic oxygen between the +2 and +1 sugars |
| 4GA500 | whole residue (21 atoms) | +1 sugar (leaving group), its O4 is the glycosidic oxygen of the scissile bond |
| 4GA501 | whole residue (21 atoms) | -1 sugar, its C1 is attacked by Asp196 |
| WAT511 | whole molecule | water bridging Glu232 and Asp299 in the active site |

<br/>

- <p align="justify">Five bonds are cut: C&alpha;-C&beta; of Asp196 and Asp299, C&beta;-C&gamma; of Glu232, and two glycosidic C-O bonds on both sides of the QM sugars, C4(499)-O4(499) and O4(501)-C1(502). The C-O cuts are polar, which the guidelines advise against, but they are a common compromise for oligosaccharides: the alternative is to include the neighbouring glucose units whole. Both are at least four bonds away from the reacting atoms (C1 of 501 and O4 of 500), and each gives a capping hydrogen on an oxygen atom, which turns the boundary oxygen into a hydroxyl group.</p>

- <p align="justify">The charge of the region (-2) comes from the two aspartates: Glu232 is neutral and the sugars and the water have no net charge. This is the &DFT CHARGE of the QM/MM inputs. With GTH pseudopotentials the region and its five capping hydrogens have 216 valence electrons, an even number, consistent with the closed-shell singlet (MULTIPLICITY 1) of the inputs.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> From a PyMOL selection to VMD and Amber </b></p> </h3>

<br/>

<p align="justify">Picking atoms by clicking is easiest in PyMOL. Write a PDB file from the topology and the coordinates with ambpdb (so that the residues keep the sequential numbering of the topology), open it in PyMOL, select the QM atoms and save the selection as MOL2:</p>

<pre style="color: white; background-color: black;">
ambpdb -p amylase.prmtop -c frame_xray.inpcrd > amylase.pdb
pymol amylase.pdb
PyMOL> save hl_65atoms.mol2, sele
</pre>

<br/>

<p align="justify">mol2_qmsel.sh reads the MOL2 file (the residue labels, such as ASP196, are RESNAME+RESID) and writes the same selection in the syntax of VMD, PyMOL and Amber masks, residue by residue. With -o the selections are also saved to files:</p>

<pre style="color: white; background-color: black;">
./mol2_qmsel.sh -o qm_sel hl_65atoms.mol2
</pre>

<pre style="color: white; background-color: black;">
RESIDUE LIST (65 atoms):
ASP196
GLH232
ASP299
4GA499
4GA500
4GA501
WAT511

VMD selection:
(name CB CG OD1 OD2 HB2 HB3 and resname ASP and resid 196) or (name CG CD OE1 OE2 HG2 HG3 HE2 and resname GLH and resid 232) or (name CB CG OD1 OD2 HB2 HB3 and resname ASP and resid 299) or (name O4 and resname 4GA and resid 499) or (name C1 C2 O2 C3 O3 C4 O4 C5 O5 C6 O6 H1 H2 H3 H4 H5 H2O H3O H61 H62 H6O and resname 4GA and resid 500) or (name C1 C2 O2 C3 O3 C4 O4 C5 O5 C6 O6 H1 H2 H3 H4 H5 H2O H3O H61 H62 H6O and resname 4GA and resid 501) or (name O H1 H2 and resname WAT and resid 511)

PYMOL selection:
(name CB+CG+OD1+OD2+HB2+HB3 &amp; resn ASP &amp; resi 196) | (name CG+CD+OE1+OE2+HG2+HG3+HE2 &amp; resn GLH &amp; resi 232) | ...

AMBER MASK selection:
(@CB,CG,OD1,OD2,HB2,HB3&amp;:ASP&amp;:196)|(@CG,CD,OE1,OE2,HG2,HG3,HE2&amp;:GLH&amp;:232)|(@CB,CG,OD1,OD2,HB2,HB3&amp;:ASP&amp;:299)|...

Selections written to qm_sel.vmd, qm_sel.pml and qm_sel.mask
</pre>

<br/>

- <p align="justify">The VMD selection (qm_sel.vmd) is the input of the next section. The Amber mask works in ParmEd and cpptraj, e.g. to list the atoms (printDetails) or to sum their MM charges, which is not an integer, as discussed in Section IV:</p>

<pre style="color: white; background-color: black;">
echo "netCharge $(cat qm_sel.mask)" | parmed -n amylase.prmtop
</pre>

<pre style="color: white; background-color: black;">
> The net charge of (@CB,CG,OD1,OD2,HB2,HB3&amp;:ASP&amp;:196)|...|(@O,H1,H2&amp;:WAT&amp;:511) is -2.2063
</pre>

<br/>

- <p align="justify">Whatever the tool, write the final selection down together with the topology it applies to: the atom numbers in the CP2K inputs are only meaningful for that topology.</p>

<br/>
