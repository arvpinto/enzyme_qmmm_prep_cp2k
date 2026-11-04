<br>
<h2> <p align="center"> <b>V - MM minimization in CP2K</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
amylase.prmtop + frame_xray.rst7 + forceeval_subsys.cp2k.inc  ->  mm_em.inp (FIST, GEO_OPT)  ->  MM_EM-1.restart
</pre>

---

<br/>

<p align="justify">The first CP2K run is a classical minimization of the whole system with the force field (METHOD FIST). It relaxes the contacts left by tleap and by the new box, and it is the first test that CP2K reads the topology correctly. Its input and outputs are in <a href="https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/tree/main/mm_em" target="_blank">mm_em</a>.</p>

<br/>

<h3> <p align="center"><b> The input </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
@SET PARM amylase.prmtop
@SET COORD frame_xray.rst7
&GLOBAL
  PROJECT MM_EM
  PRINT_LEVEL LOW
  RUN_TYPE GEO_OPT
&END GLOBAL

&FORCE_EVAL
  METHOD FIST
  STRESS_TENSOR ANALYTICAL
  &MM
    &FORCEFIELD
      DO_NONBONDED
      EI_SCALE14 0.8333333
      IGNORE_MISSING_CRITICAL_PARAMS FALSE
      PARMTYPE AMBER
      PARM_FILE_NAME $PARM
      SHIFT_CUTOFF TRUE
      VDW_SCALE14 0.5
      &SPLINE
        EMAX_SPLINE 1.0E8
        RCUT_NB [angstrom] 10
      &END SPLINE
    &END FORCEFIELD
    ...
    &POISSON
      PERIODIC XYZ
      POISSON_SOLVER PERIODIC
      &EWALD
        ALPHA .40
        EWALD_TYPE SPME
        GMAX 80
        O_SPLINE 6
        RCUT 10.0
        ...
      &END EWALD
    &END POISSON
  &END MM
  &RESCALE_FORCES
    MAX_FORCE 4.21648784E-2
  &END RESCALE_FORCES
  &SUBSYS
    @INCLUDE forceeval_subsys.cp2k.inc
    &TOPOLOGY
      CONN_FILE_FORMAT AMBER
      CONN_FILE_NAME $PARM
      COORD_FILE_FORMAT CRD
      COORD_FILE_NAME $COORD
      ...
    &END TOPOLOGY
    &KIND NA+
     ELEMENT Na
    &END KIND
    &KIND CL-
     ELEMENT Cl
    &END KIND
  &END SUBSYS
&END FORCE_EVAL

&MOTION
  &GEO_OPT
    MAX_DR    1.0E-02
    MAX_FORCE 1.0E-02
    MAX_ITER 1000
    OPTIMIZER LBFGS
    RMS_DR    5.0E-03
    RMS_FORCE 5.0E-03
    TYPE MINIMIZATION
  &END GEO_OPT
  ...
&END MOTION
</pre>

<br/>

- <p align="justify">The force field, the connectivity and the charges come from amylase.prmtop (PARMTYPE AMBER and CONN_FILE_FORMAT AMBER), the coordinates from frame_xray.rst7, read as an Amber CRD file, and the cell from forceeval_subsys.cp2k.inc. The cell has to be given: CP2K reads the box of a CRD file only when it also has velocities, and the rst7 written by tleap has none. The &KIND sections give the element of the two ions whose kind names CP2K cannot map (Section I).</p>

- <p align="justify">EI_SCALE14 and VDW_SCALE14 are the inverse of the Amber SCEE (1.2) and SCNB (2.0) factors of ff14SB. CP2K applies them to every 1-4 pair and does not read the per-dihedral SCEE_SCALE_FACTOR and SCNB_SCALE_FACTOR sections of the prmtop. GLYCAM uses no 1-4 scaling (SCEE = SCNB = 1.0, as stored in the topology for its dihedrals), so in CP2K the 1-4 interactions of the oligosaccharide are scaled as in the protein. This affects only the MM sugars (the QM sugars are described by DFT), but it changes the MM conformational energies of the carbohydrate, and it is a reason why CP2K and Amber energies of this system will not agree exactly.</p>

- <p align="justify">SHIFT_CUTOFF TRUE shifts every real-space pair interaction to zero at the cutoff (10 &Aring;), which improves energy conservation in MD but adds a constant per pair: CP2K MM energies are then not directly comparable with Amber. Keep the same setting in all the calculations you compare.</p>

- <p align="justify">The electrostatics use SPME with ALPHA 0.40, a real-space cutoff of 10 &Aring; and GMAX 80. GMAX 80 is below the rule of thumb of one grid point per &Aring;ngstr&ouml;m for the 103 &Aring; edge of this box; if you validate the MM energy against Amber, test a finer grid (e.g. GMAX 91 103 91).</p>

- <p align="justify">RESCALE_FORCES caps the force on any atom at 0.042 hartree/bohr (about 50 kcal mol<sup>-1</sup> &Aring;<sup>-1</sup>), so that a bad contact in the starting structure does not throw atoms across the box in the first steps.</p>

- <p align="justify">The convergence criteria are loose (forces of 10<sup>-2</sup> hartree/bohr): the goal is to remove the strain, not to reach a minimum. A frame is written to MM_EM-pos-1.dcd every 100 steps, and the restart file MM_EM-1.restart, the starting point of the QM/MM minimization, every 500 steps and at the end.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Running the minimization </b></p> </h3>

<br/>

<p align="justify">Each CP2K run reads files from system_prep and from the directory of the previous run. Copy them into the run directory (preferably a copy of it, in your scratch space) and run CP2K there:</p>

<pre style="color: white; background-color: black;">
cp -r ~/enzyme_qmmm_prep_cp2k/mm_em ~/scratch/ && cd ~/scratch/mm_em
cp ~/enzyme_qmmm_prep_cp2k/system_prep/{amylase.prmtop,frame_xray.rst7,forceeval_subsys.cp2k.inc} .
mpirun -np 64 cp2k.psmp -i mm_em.inp -o mm_em.out
</pre>

<p align="justify">or, from the repository directory, let run_cp2k.sh copy the files and run the stages in order, in a work directory that mirrors the repository:</p>

<pre style="color: white; background-color: black;">
./scripts/run_cp2k.sh -h
./scripts/run_cp2k.sh -w ~/scratch/amylase -s mm_em -n 64
./scripts/run_cp2k.sh -w ~/scratch/amylase -n 128 -l srun      # all stages, inside a Slurm job
./scripts/run_cp2k.sh -w ~/scratch/amylase -d                  # print what would be done
</pre>

<br/>

<p align="justify">In the example, the minimization ran the 1000 steps of MAX_ITER in 4 minutes on 64 MPI ranks, and stopped with "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED". With criteria this loose, and the energy down from -289.5 to -433.8 hartree, that is enough to start the QM/MM minimization; the cp2k_opt function of cp2k_funcs prints the convergence criteria at every step:</p>

<pre style="color: white; background-color: black;">
source cp2k_funcs
cp2k_opt mm_em.out | tail -6
</pre>

<br/>

- <p align="justify">Read the beginning of mm_em.out: the number of atoms, kinds and molecules, the cell and the charge of the system (CHARGE_INFO| Total Charge of the Classical System: -0.000000). Three warnings are expected:</p>

<pre style="color: white; background-color: black;">
 *** WARNING in topology_amber.F:235 :: No VELOCITY information found in ***
 *** CRD file. Ignoring BOX information. Please provide the BOX          ***
 *** information directly from the main CP2K input!                      ***
 *** WARNING in topology_amber.F:260 :: BOX information missing in CRD ***
 FORCEFIELD| WARNING: A non Critical ForceField parameter is missing! CP2K GOES!
 FORCEFIELD| Non critical parameters are:Urey-Bradley,Torsions,Impropers, Opbends and 1-4
</pre>

<p align="justify">The first two are the reason for forceeval_subsys.cp2k.inc. The last one is expected because Amber has no Urey-Bradley terms; activate &MM/&PRINT/&FF_INFO to list the missing terms and make sure that no bond, angle or torsion is among them.</p>

- <p align="justify">Do not write an inline comment after PARM_FILE_NAME or CONN_FILE_NAME: CP2K takes it as part of the file name and stops with "Force field file missing".</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Looking at the trajectory </b></p> </h3>

<br/>

<div class="qm-movie">
  <video src="assets/movies/mm_em_qm.mp4" autoplay loop muted playsinline controls preload="metadata" aria-label="QM region during the MM minimization"></video>
</div>

<p align="center"><i>The QM region during the MM minimization, from the starting structure (frame_xray.pdb) to the minimized one (MM_EM-1.restart). Yellow dashes: Asp196 OD1-C1 and Glu232 HE2-O4. The run writes few frames, so they are joined by cubic interpolation and the motion between them is only illustrative.</i></p>

<br/>

<p align="justify">The trajectories of the example are attached to the <a href="https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/releases/latest" target="_blank">latest release</a> of the repository. Download MM_EM-pos-1.dcd into mm_em and load it with the topology:</p>

<pre style="color: white; background-color: black;">
cd mm_em
wget https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/releases/latest/download/MM_EM-pos-1.dcd
vmd -parm7 ../system_prep/amylase.prmtop -dcd MM_EM-pos-1.dcd
</pre>

<p align="justify">Use the selections of Section I (all, protein, resname 0GA 4GA ROH) and check that the protein, the oligosaccharide and the ions stay in place, and that the water at the faces of the box closes the gaps left by the new box. For the active site, vmd_cp2k-repr.tcl draws the QM region with dynamic bonds (Section VIII):</p>

<pre style="color: white; background-color: black;">
vmd -e ../scripts/vmd_cp2k-repr.tcl -args ../system_prep/amylase.prmtop MM_EM-pos-1.dcd \
    -qminc ../system_prep/forceeval_qmmm.cp2k.inc -cartoon
</pre>

<br/>
