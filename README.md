<br>
<p align="justify"><b>This tutorial prepares, step by step, a QM/MM model of an enzyme-substrate complex for CP2K and takes it all the way to a reactive simulation: starting from a PDB frame, it builds the Amber topology with tleap, makes it readable by CP2K, defines the QM region and the &QMMM sections, corrects the MM charges for electrostatic embedding, minimizes and equilibrates the system with QM/MM, and ends with a steered QM/MM MD run in which the glycosidic bond of the substrate breaks while a proton is transferred from the catalytic acid.</b></p>

<p align="justify">Most of the information here came from the <a href="https://www.cp2k.org/" target="_blank">CP2K site</a>, the <a href="https://manual.cp2k.org/" target="_blank">CP2K manual</a>, the <a href="https://ambermd.org/tutorials/" target="_blank">Amber tutorials</a>, the <a href="https://groups.google.com/g/cp2k" target="_blank">CP2K discussion group</a> and the BioExcel QM/MM material listed at the bottom of this page. You should also follow the good practices found in the literature for the specific chemistry of your system.</p>

<br>

<div align="center">
    <img src="workflow.png" width="90%">
</div>

<br>

---

<br>

<h3><p align="center"><b> Example system </b></p></h3>

<p align="justify">All the commands and files of the tutorial come from a worked example: human pancreatic &alpha;-amylase (HPA), a retaining glycoside hydrolase of family GH13, with a maltopentaose chain bound in the active site. In retaining &alpha;-amylases the first step of the reaction (glycosylation) breaks the &alpha;-1,4 glycosidic bond between the sugar in subsite -1 and the one in subsite +1: the general acid/base (Glu233 in HPA) protonates the glycosidic oxygen, the nucleophile (Asp197) attacks the anomeric carbon of the -1 sugar and a third carboxylate (Asp300) binds the hydroxyl groups of that sugar. The structure starts at Tyr2, so these residues are Asp196, Glu232 and Asp299 in the topology, and Glu232 is protonated (GLH) to act as the acid.</p>

<p align="justify">The model is a frame of a classical MD simulation of the solvated complex (frame_xray.pdb): 67 992 atoms in an orthorhombic box of 88.6 x 102.3 x 90.2 &Aring;, with 20 076 TIP3P water molecules, one Ca<sup>2+</sup>, one Cl<sup>-</sup> and one Na<sup>+</sup> ion. The protein is described with ff99SB and the oligosaccharide with GLYCAM_06h (residues ROH, 4GA and 0GA, numbered 498 to 503 in the topology). The QM region has 65 atoms (PBE/DZVP, GPW) and the steered MD drives the collective variable CV = d(O4-C1) - d(H<sub>Glu232</sub>-O4), which goes from negative (glycosidic bond intact, proton on Glu232) to positive (bond broken, proton on the leaving sugar).</p>

<p align="justify">The steps are general and apply to any enzyme-substrate complex parameterized with Amber-type force fields. Replace the residue numbers, atom numbers and file names with those of your system.</p>

<br>

<h3><p align="center"><b> Requirements </b></p></h3>

- <p align="justify"><a href="https://ambermd.org/AmberTools.php" target="_blank">AmberTools</a> (tleap and ParmEd).</p>
- <p align="justify"><a href="https://www.cp2k.org/" target="_blank">CP2K</a> compiled with MPI (cp2k.psmp), with its data directory (BASIS_SET, POTENTIAL and MM_POTENTIAL files) reachable through CP2K_DATA_DIR. The QM/MM steps need a cluster: every MD step is a DFT calculation of the QM region in the field of about 68 000 MM atoms.</p>
- <p align="justify"><a href="https://www.ks.uiuc.edu/Research/vmd/" target="_blank">VMD</a> to build the &QMMM include files and to look at the trajectories, and optionally <a href="https://pymol.org/" target="_blank">PyMOL</a> to pick the QM region.</p>
- <p align="justify">Python 3 with ParmEd, NumPy and matplotlib (pip install parmed numpy matplotlib) for the Python scripts, and bash for the shell scripts.</p>

<br>

<h3><p align="center"><b> Files </b></p></h3>

| File | Description |
| :--- | :--- |
| <a href="scripts/prmtop_fix_types.py">prmtop_fix_types.py</a> | Checks the atom types of a prmtop against the way CP2K reads them, and renames those that would be merged or misread |
| <a href="scripts/mol2_qmsel.sh">mol2_qmsel.sh</a> | Converts a MOL2 file of the QM region (e.g. saved from PyMOL) into VMD, PyMOL and Amber mask selections |
| <a href="scripts/vmd_cp2k-qmmm.tcl">vmd_cp2k-qmmm.tcl</a> | VMD script that writes the &CELL, &QM_KIND and &LINK include files and the ParmEd script of the charge correction |
| <a href="scripts/run_cp2k.sh">run_cp2k.sh</a> | Runs the CP2K stages in order, checking the files that each one needs |
| <a href="scripts/cp2k_funcs">cp2k_funcs</a> | Bash functions to read CP2K outputs (optimization convergence, MD energies, CV, QM atom serials) |
| <a href="scripts/vmd_cp2k-repr.tcl">vmd_cp2k-repr.tcl</a> | VMD script that loads CP2K trajectories and draws the QM region with dynamic bonds |
| <a href="scripts/plot_smd.py">plot_smd.py</a> | Plots the CV, the restraint target, the work and the potential energy of the steered MD run |
| <a href="AMBER_INPUTS.zip">AMBER_INPUTS.zip</a> | frame_xray.pdb, leap_amyl.in and the MOL2 file of the QM region |
| <a href="CP2K_INPUTS.zip">CP2K_INPUTS.zip</a> | CP2K inputs (em, em-me_gpw, em-ee_gpw, react_nvt-ee_gpw, smd-nvt-ee_gpw) and the files generated for the example |

<br>

<p align="center"><b> Once you have a PDB structure of your enzyme-substrate complex with Amber residue names, you are ready to start this tutorial. </b></p>

---

### Useful Links:

- <a href="https://manual.cp2k.org/trunk/methods/qm_mm/builtin.html" target="_blank">CP2K manual: QM/MM with the built-in QM/MM driver</a>
- <a href="https://docs.bioexcel.eu/2020_06_09_online_ambertools4cp2k/" target="_blank">Preparing to run biomolecular QM/MM simulations with CP2K using AmberTools</a>
- <a href="https://github.com/bioexcel/cp2K_qmmm_tutorials_for_biological_simulations" target="_blank">QM/MM tutorials adapted to CP2K</a>
- <a href="https://docs.bioexcel.eu/qmmm_simulation_bpg/en/main/index.html" target="_blank">Best Practices in QM/MM Simulation of Biomolecular Systems</a>
- <a href="https://glycam.org/" target="_blank">GLYCAM force field for carbohydrates</a>
- <a href="https://arvpinto.github.io/enzyme_neb_cp2k/" target="_blank">Tutorial for TS searches of enzymatic reaction mechanisms using QM/MM methods in CP2K</a>
