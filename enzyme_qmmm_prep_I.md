<br>
<h2> <p align="center"> <b>I - System setup with tleap and atom types for CP2K</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
frame_xray.pdb  ->  tleap (ff14SB + GLYCAM_06j-1 + TIP3P)  ->  amylase_prep.prmtop + frame_xray.rst7  ->  prmtop_fix_types.py  ->  amylase.prmtop
</pre>

---

<br/>

<p align="justify">In this section we look at the starting structure, build the Amber topology (prmtop) and coordinates (rst7) of the solvated complex, and check that CP2K will read the topology as Amber does. The second step is short, but skipping it can give wrong MM electrostatics without any error message.</p>

<p align="justify">The files of Sections I to IV, inputs and outputs, are in <a href="https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/tree/main/system_prep" target="_blank">system_prep</a>. To follow along, work in a new directory with the three inputs (frame_xray.pdb, leap_amyl.in and hl_65atoms.mol2) and the scripts, and compare what you get with the files of system_prep. The same directory is used for the CP2K runs of Sections V to VIII:</p>

<pre style="color: white; background-color: black;">
git clone https://github.com/arvpinto/enzyme_qmmm_prep_cp2k.git
mkdir amylase_setup && cd amylase_setup
cp ../enzyme_qmmm_prep_cp2k/system_prep/{frame_xray.pdb,leap_amyl.in,hl_65atoms.mol2} .
cp ../enzyme_qmmm_prep_cp2k/scripts/* .
</pre>

<br/>

<h3> <p align="center"><b> The starting structure </b></p> </h3>

<br/>

<p align="justify">frame_xray.pdb is a frame of a classical MD simulation of the complex, so it already has everything tleap needs: hydrogen atoms with Amber names, water molecules and ions, and residue names that encode the protonation states and the carbohydrate. Open it in VMD:</p>

<pre style="color: white; background-color: black;">
vmd frame_xray.pdb
</pre>

<p align="justify">and, in Graphics > Representations, look at it with different selections:</p>

| Selection | What to check |
| :--- | :--- |
| all | the water box around the protein, with the ions (67 992 atoms) |
| protein | 495 residues, Tyr2 to Leu496 in HPA numbering (1 to 495 here) |
| resname 0GA 4GA ROH | the maltopentaose in the active site cleft: ROH 498 caps the reducing end, 4GA 499 to 502 are &alpha;-D-glucopyranoses linked at O4 to the next unit, 0GA 503 is the non-reducing end |
| resname GLH or resid 196 299 | the catalytic residues around the -1 and +1 sugars (4GA501 and 4GA500) |

<br/>

- <p align="justify">The residue names define the protonation states: GLH232 (Glu233 in HPA numbering) is neutral and protonated at OE2, because it is the general acid of the glycosylation step; HIE are histidines protonated at N&epsilon;2; CYX are the cysteines of the five disulfide bonds (27-85, 69-114, 140-159, 377-383, 449-461). CA2, Cl- and Na+ are the Ca<sup>2+</sup>, Cl<sup>-</sup> and Na<sup>+</sup> ions.</p>

- <p align="justify">If you start from a crystal structure instead, the usual preparation applies before this point: missing atoms and loops, protonation states (PROPKA, H++ and inspection of the hydrogen bond network, in particular of the catalytic residues), the GLYCAM residue names of the carbohydrate (the <a href="https://glycam.org/" target="_blank">GLYCAM web tools</a> build and name oligosaccharides), solvation, neutralization and a classical MD equilibration. The protonation state of the catalytic residues has to match the mechanism you want to study: here the acid/base Glu232 must carry the proton that it will transfer, and the nucleophile Asp196 must be deprotonated.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Building the topology with tleap </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
# leap_amyl.in - Amber topology and coordinates of the alpha-amylase/oligosaccharide model
# Usage: tleap -f leap_amyl.in > leap_amyl.log
source leaprc.protein.ff14SB
source leaprc.GLYCAM_06j-1
source leaprc.water.tip3p

amyl = loadpdb frame_xray.pdb

bond amyl.27.SG amyl.85.SG
bond amyl.69.SG amyl.114.SG
bond amyl.140.SG amyl.159.SG
bond amyl.377.SG amyl.383.SG
bond amyl.449.SG amyl.461.SG

alignaxes amyl

setBox amyl vdw {0.0 0.0 0.0}

# the atom types are checked/renamed for CP2K afterwards (prmtop_fix_types.py -> amylase.prmtop)
saveamberparm amyl amylase_prep.prmtop frame_xray.rst7

quit
</pre>

<pre style="color: white; background-color: black;">
tleap -f leap_amyl.in
</pre>

<br/>

- <p align="justify">The protein is described with ff14SB, the carbohydrate with GLYCAM_06j-1 and the water with TIP3P; leaprc.water.tip3p also loads the ion parameters for TIP3P (Joung-Cheatham for Na<sup>+</sup> and Cl<sup>-</sup>, Li/Merz 12-6 for Ca<sup>2+</sup>). The disulfide bonds are defined with the bond command; do not also keep CONECT records for them in the PDB file, or tleap will add each bond twice.</p>

- <p align="justify">tleap writes its log to leap.log. The end of the log should read:</p>

<pre style="color: white; background-color: black;">
Box dimensions:  90.757909 103.275752 91.248328
...
Exiting LEaP: Errors = 0; Warnings = 0; Notes = 0.
</pre>

<p align="justify">The "No sp2 improper torsion term for H-H-N3-H" and "C8-H-N3-H" messages printed while the topology is built (92 here) refer to the NH<sub>3</sub><sup>+</sup> groups of the N-terminus and of the lysines, which have no improper torsion in ff14SB: they are expected and not counted as warnings. Look instead for missing parameters, atoms added that should have been in the PDB, unexpected bonds and the total charge (0 here).</p>

- <p align="justify">alignaxes rotates the system so that its principal axes lie along x, y and z, and setBox vdw sets an orthorhombic box that encloses all atoms (90.758 x 103.276 x 91.248 &Aring;, the last line of frame_xray.rst7). This is convenient for a frame without box information, but if the frame comes from a periodic simulation the original box is a better choice: rotating the frame or changing the box breaks the continuity of the water across the periodic boundary, and the first minimization has to repair it. In that case remove alignaxes and set the box of the frame with set amyl box {a b c}.</p>

- <p align="justify">Open the topology and the coordinates in VMD and repeat the checks above, now with the topology (bonds, disulfides, the position of the oligosaccharide after alignaxes):</p>

<pre style="color: white; background-color: black;">
vmd -parm7 amylase_prep.prmtop -rst7 frame_xray.rst7
</pre>

<br/>

- <p align="justify">From now on, take every residue and atom number from the topology, never from an original crystal structure. In VMD, resid follows the sequential residue numbering of the prmtop, serial is the 1-based atom number used by CP2K, and index is serial - 1.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Atom types that CP2K can read </b></p> </h3>

<br/>

<p align="justify">CP2K reads the Amber topology directly (PARMTYPE AMBER and CONN_FILE_FORMAT AMBER), but it does not use the atom types as Amber does. It creates one atomic kind per atom type and charge, and only the names of these kinds tell the atoms apart:</p>

- <p align="justify">the type names are converted to upper case, so GLYCAM Ho/Oh (sugar hydroxyl) and ff14SB HO/OH (Ser, Thr, Tyr and GLH hydroxyl) fall into the same group of kinds;</p>

- <p align="justify">when a type has several charges, its kinds are named by appending a counter: the amide hydrogens (type H, on the backbone and on the Asn, Gln and Arg side chains, with 12 different charges here) become H1, H2, ..., H12. The topology also has types called H1, H2 and H5 (aliphatic H next to an electronegative atom in the protein and in GLYCAM, the anomeric H of GLYCAM and the H on C&epsilon;1 of His), whose kinds are H11, H12, ..., H131, H2 and H5. Atoms of two different types then share one kind, with a single charge. Released versions of CP2K, including the 2023.1 version used for this example and the latest release (2026.2), do not stop or warn: the electrostatics are simply wrong. This was fixed in the development version of CP2K in October 2026 (<a href="https://github.com/cp2k/cp2k/pull/6191" target="_blank">cp2k/cp2k#6191</a>);</p>

- <p align="justify">the element of each kind is guessed from its name, with the trailing digits removed, using built-in lists of Amber, CHARMM and GROMOS types and the periodic table. The calcium type, C0, becomes carbon, and Na+ and Cl- are not recognized at all (the inputs of this tutorial give them a &KIND section with their ELEMENT).</p>

<br/>

<p align="justify">prmtop_fix_types.py reproduces the kind names that CP2K builds, reports these problems and writes a topology in which the offending types are renamed. The new names are taken from the type lists that CP2K already knows, so that the element is guessed correctly and no extra &KIND sections are needed. Only the labels in AMBER_ATOM_TYPE change: the parameters are indexed by number in the prmtop, so Amber gives exactly the same energy with both files.</p>

<pre style="color: white; background-color: black;">
python3 prmtop_fix_types.py -p amylase_prep.prmtop -o amylase.prmtop
</pre>

<pre style="color: white; background-color: black;">
prmtop_fix_types.py v2.0
amylase_prep.prmtop: 43 atom types, 259 CP2K kinds
  type   element    atoms  charges  CP2K kinds
  N3     N             23        2  N31 ... N32
  H      H            818       12  H1 ... H12
  ...
  H1     H            828       31  H11 ... H131
  ...
  C8     C            172        7  C81 ... C87
  ...
  H5     H             11        1  H5
  ...
  C0     Ca             1        1  C0
  Cl-    Cl             1        1  CL-
  Ho     H             17        5  HO1 ... HO10
  Oh     O             17        5  OH1 ... OH10
  Cg     C             30        6  CG1 ... CG6
  H2     H              5        1  H2
  ...
  kind names shared by different types: 4 (H2 &lt;- H/H2, H5 &lt;- H/H5, H11 &lt;- H/H1, H12 &lt;- H/H1)

Types whose element CP2K cannot guess (give them a &KIND section with ELEMENT in &SUBSYS):
  Cl-    element Cl  kind name(s): CL-
  Na+    element Na  kind name(s): NA+

old    new   element    atoms  reason
N3     NT    N             23  name can be produced by the kind counter of N
H1     HZ    H            828  name can be produced by the kind counter of H
C8     CZ    C            172  name can be produced by the kind counter of C
N2     NL    N             84  name can be produced by the kind counter of N
H5     HN    H             11  name can be produced by the kind counter of H
H4     HL    H             28  name can be produced by the kind counter of H
O2     OA    O            106  name can be produced by the kind counter of O
H2     HR    H              5  name can be produced by the kind counter of H
Oh     OD    O             17  same upper-case name as OH
Ho     HT    H             17  same upper-case name as HO
C0     CAL   Ca             1  CP2K reads it as C, element is Ca

amylase.prmtop: 43 atom types, 259 CP2K kinds
  kind names shared by different types: none

11 type(s) renamed. Topology written to amylase.prmtop; table in type_renaming.dat
</pre>

<br/>

- <p align="justify">Four kind names are shared in the tleap topology: the GLYCAM anomeric H (H2) and the His H5 each share a kind with amide hydrogens, and two kinds of the aliphatic H1 hydrogens do the same. The script renames every type whose name can be produced by the counter of another type (N2, N3, O2, C8, H1, H2, H4, H5), even when there is no clash yet: the charge correction of Section IV creates new charges, and with them new counters.</p>

- <p align="justify">Use --dry-run to check a topology without writing anything (Section IV does it for the corrected topology), --map OLD:NEW to choose the new names yourself and --keep to protect a type. The table of renamed types is saved to type_renaming.dat; keep it with your inputs, since the new type names appear in the CP2K output and in VMD.</p>

- <p align="justify">From now on, amylase.prmtop is the topology of every MM and mechanical embedding calculation, together with frame_xray.rst7.</p>

<br/>
