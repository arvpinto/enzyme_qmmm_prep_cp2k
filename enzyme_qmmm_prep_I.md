<br>
<h2> <p align="center"> <b>I - System setup with tleap and atom types for CP2K</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
frame_xray.pdb  ->  tleap (ff99SB + GLYCAM_06h)  ->  amylase_leap.prmtop + frame_xray.inpcrd  ->  prmtop_fix_types.py  ->  amylase.prmtop
</pre>

---

<br/>

<p align="justify">In this section we build the Amber topology (prmtop) and coordinates (inpcrd) of the solvated complex, and then check that CP2K will read the topology as Amber does. The second step is short, but skipping it can give wrong MM electrostatics without any error message.</p>

<br/>

<h3> <p align="center"><b> The starting structure </b></p> </h3>

<br/>

- <p align="justify">frame_xray.pdb (in <a href="AMBER_INPUTS.zip">AMBER_INPUTS.zip</a>) is a frame of a classical MD simulation of the complex, so it already has everything tleap needs: hydrogen atoms with Amber names, water molecules and ions, and residue names that encode the protonation states and the carbohydrate:</p>

| Residue names | Meaning |
| :--- | :--- |
| GLH232 | Glu232 (Glu233 in HPA numbering) neutral, protonated at OE2: the general acid of the glycosylation step |
| HIE | histidines protonated at N&epsilon;2 |
| CYX | cysteines in disulfide bonds (27-85, 69-114, 140-159, 377-383, 449-461) |
| ROH 498, 4GA 499-502, 0GA 503 | GLYCAM maltopentaose: ROH caps the reducing end, 4GA is an &alpha;-D-glucopyranose linked at O4 to the next unit, 0GA the non-reducing end |
| CA2, Cl-, Na+, WAT | Ca<sup>2+</sup>, Cl<sup>-</sup>, Na<sup>+</sup> and TIP3P water |

<br/>

- <p align="justify">If you start from a crystal structure instead, the usual preparation applies before this point: missing atoms and loops, protonation states (PROPKA, H++ and inspection of the hydrogen bond network, in particular of the catalytic residues), the GLYCAM residue names of the carbohydrate (the <a href="https://glycam.org/" target="_blank">GLYCAM web tools</a> build and name oligosaccharides), solvation, neutralization and a classical MD equilibration. The protonation state of the catalytic residues has to match the mechanism you want to study: here the acid/base Glu232 must carry the proton that it will transfer, and the nucleophile Asp196 must be deprotonated.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Building the topology with tleap </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
# leap_amyl.in
source oldff/leaprc.ff99SB
source oldff/leaprc.GLYCAM_06h

loadamberparams frcmod.ionsff99_tip3p

amyl = loadpdb frame_xray.pdb

bond amyl.27.SG amyl.85.SG
bond amyl.69.SG amyl.114.SG
bond amyl.140.SG amyl.159.SG
bond amyl.377.SG amyl.383.SG
bond amyl.449.SG amyl.461.SG

alignaxes amyl

setBox amyl vdw {0.0 0.0 0.0}

saveamberparm amyl amylase_leap.prmtop frame_xray.inpcrd

quit
</pre>

<pre style="color: white; background-color: black;">
tleap -f leap_amyl.in > leap_amyl.log
</pre>

<br/>

- <p align="justify">The protein is described with ff99SB and the carbohydrate with GLYCAM_06h, with the TIP3P water model and its ion parameters. These are the force fields of the MD simulation that produced the frame: always keep the force field of the equilibration. For a new system, ff14SB or ff19SB and GLYCAM_06j are the current choices.</p>

- <p align="justify">The disulfide bonds are defined with the bond command (the residues are already called CYX in the PDB). Do not also keep CONECT records for them in the PDB file, or tleap will add each bond twice.</p>

- <p align="justify">alignaxes rotates the system so that its principal axes lie along x, y and z, and setBox vdw sets an orthorhombic box that encloses all atoms (here 88.589 x 102.293 x 90.200 &Aring;). This is convenient for a frame without box information, but if the frame comes from a periodic simulation the original box is a better choice: rotating the frame or changing the box breaks the continuity of the water across the periodic boundary, and the first minimization has to repair it. In that case remove alignaxes and set the box of the frame with set amyl box {a b c}.</p>

- <p align="justify">Read the whole leap_amyl.log: missing parameters, atoms added that should have been in the PDB, unexpected bonds and the total charge (0 here). Then open the topology and the coordinates in VMD and check the protonation states, the disulfides and the position of the oligosaccharide.</p>

<pre style="color: white; background-color: black;">
vmd -parm7 amylase_leap.prmtop -rst7 frame_xray.inpcrd
</pre>

<br/>

- <p align="justify">From now on, take every residue and atom number from the topology (VMD resid and serial, or ParmEd), never from an original crystal structure. In VMD, resid follows the sequential residue numbering of the prmtop, serial is the 1-based atom number used by CP2K, and index is serial - 1.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Atom types that CP2K can read </b></p> </h3>

<br/>

<p align="justify">CP2K reads the Amber topology directly (PARMTYPE AMBER and CONN_FILE_FORMAT AMBER), but it does not use the atom types as Amber does. It creates one atomic kind per atom type and charge, and only the names of these kinds tell the atoms apart:</p>

- <p align="justify">the type names are converted to upper case, so GLYCAM Ho/Oh (sugar hydroxyl) and ff99SB HO/OH (Ser, Thr, Tyr and GLH hydroxyl) fall into the same group of kinds;</p>

- <p align="justify">when a type has several charges, its kinds are named by appending a counter: the amide hydrogens (type H, on the backbone and on the Asn, Gln and Arg side chains, with 12 different charges here) become H1, H2, ..., H12. The topology also has types called H1, H2 and H5 (aliphatic H next to an electronegative atom in the protein and in GLYCAM, the anomeric H of GLYCAM and the H on C&epsilon;1 of His), whose kinds are H11, H12, ..., H131, H2 and H5. Atoms of two different types then share one kind, with a single charge. Released versions of CP2K (up to 2026.2) do not stop or warn: the electrostatics are simply wrong. This was fixed in the development version of CP2K in October 2026 (<a href="https://github.com/cp2k/cp2k/pull/6191" target="_blank">cp2k/cp2k#6191</a>);</p>

- <p align="justify">the element of each kind is guessed from its name, with the trailing digits removed, using built-in lists of Amber, CHARMM and GROMOS types and the periodic table. The calcium type of the old ion parameters, C0, becomes carbon, and Na+ and Cl- are not recognized at all (the inputs of this tutorial give them a &KIND section with their ELEMENT).</p>

<br/>

<p align="justify">prmtop_fix_types.py reproduces the kind names that CP2K builds, reports these problems and writes a topology in which the offending types are renamed. The new names are taken from the type lists that CP2K already knows, so that the element is guessed correctly and no extra &KIND sections are needed. Only the labels in AMBER_ATOM_TYPE change: the parameters are indexed by number in the prmtop, so Amber gives exactly the same energy with both files.</p>

<pre style="color: white; background-color: black;">
python3 prmtop_fix_types.py -p amylase_leap.prmtop -o amylase.prmtop
</pre>

<pre style="color: white; background-color: black;">
prmtop_fix_types.py v2.0
amylase_leap.prmtop: 38 atom types, 259 CP2K kinds
  type   element    atoms  charges  CP2K kinds
  N3     N             23        2  N31 ... N32
  H      H            818       12  H1 ... H12
  CT     C           1410       68  CT1 ... CT68
  ...
  H1     H            828       31  H11 ... H131
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
N2     NL    N             84  name can be produced by the kind counter of N
H5     HN    H             11  name can be produced by the kind counter of H
H4     HL    H             28  name can be produced by the kind counter of H
O2     OA    O            106  name can be produced by the kind counter of O
H2     HR    H              5  name can be produced by the kind counter of H
Oh     OD    O             17  same upper-case name as OH
Ho     HT    H             17  same upper-case name as HO
C0     CAL   Ca             1  CP2K reads it as C, element is Ca

amylase.prmtop: 38 atom types, 259 CP2K kinds
  kind names shared by different types: none

10 type(s) renamed. Topology written to amylase.prmtop; table in type_renaming.dat
</pre>

<br/>

- <p align="justify">Four kind names are shared in the original topology: the GLYCAM anomeric H (H2) and the His H5 each share a kind with backbone amide hydrogens, and two kinds of aliphatic H1 hydrogens do the same. The script renames every type whose name can be produced by the counter of another type (N2, N3, O2, H1, H2, H4, H5), even when there is no clash yet: the charge correction of Section IV creates new charges, and with them new counters. The charge-corrected topology of this tutorial would have five shared kinds with the original type names.</p>

- <p align="justify">Use --dry-run to check a topology without writing anything (for example the corrected topology of Section IV), --map OLD:NEW to choose the new names yourself and --keep to protect a type. The table of renamed types is saved to type_renaming.dat; keep it with your inputs, since the type names appear in the CP2K output and in VMD.</p>

- <p align="justify">Confirm with sander or pmemd that both topologies give the same single-point energy if you want to be sure: only the labels changed. In the rest of the tutorial, amylase.prmtop is the topology of every MM and mechanical embedding calculation.</p>

<br/>

<p align="justify">The input files of this section are in <a href="AMBER_INPUTS.zip">AMBER_INPUTS.zip</a>.</p>

<br/>
