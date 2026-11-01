<br>
<h2> <p align="center"> <b>V - MM minimization in CP2K</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
amylase.prmtop + frame_xray.inpcrd + forceeval_subsys.cp2k.inc  ->  em.inp (FIST, GEO_OPT)  ->  EM-1.restart
</pre>

---

<br/>

<p align="justify">The first CP2K run is a classical minimization of the whole system with the force field (METHOD FIST). It relaxes the contacts left by tleap and by the new box, and it is the first test that CP2K reads the topology correctly. All CP2K inputs are in <a href="CP2K_INPUTS.zip">CP2K_INPUTS.zip</a>, and the stages can be run one by one or in sequence with run_cp2k.sh.</p>

<br/>

<h3> <p align="center"><b> The input </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
@SET PARM amylase.prmtop
@SET COORD frame_xray.inpcrd
&GLOBAL
  PROJECT EM
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

- <p align="justify">The force field, the connectivity and the charges come from amylase.prmtop (PARMTYPE AMBER and CONN_FILE_FORMAT AMBER), the coordinates from the inpcrd file and the cell from forceeval_subsys.cp2k.inc. The &KIND sections give the element of the two ions whose kind names CP2K cannot map (Section I).</p>

- <p align="justify">EI_SCALE14 and VDW_SCALE14 are the inverse of the Amber SCEE (1.2) and SCNB (2.0) factors of ff99SB. CP2K applies them to every 1-4 pair and does not read the per-dihedral SCEE_SCALE_FACTOR and SCNB_SCALE_FACTOR sections of the prmtop. GLYCAM uses no 1-4 scaling (SCEE = SCNB = 1.0, as stored in the topology for its dihedrals), so in CP2K the 1-4 interactions of the oligosaccharide are scaled as in the protein. This affects only the MM sugars (the QM sugars are described by DFT), but it changes the MM conformational energies of the carbohydrate, and it is a reason why CP2K and Amber energies of this system will not agree exactly.</p>

- <p align="justify">SHIFT_CUTOFF TRUE shifts every real-space pair interaction to zero at the cutoff (10 &Aring;), which improves energy conservation in MD but adds a constant per pair: CP2K MM energies are then not directly comparable with Amber. Keep the same setting in all the calculations you compare.</p>

- <p align="justify">The electrostatics use SPME with ALPHA 0.40, a real-space cutoff of 10 &Aring; and GMAX 80. GMAX 80 is below the rule of thumb of one grid point per &Aring;ngstr&ouml;m for the 102 &Aring; edge of this box; if you validate the MM energy against Amber, test a finer grid (e.g. GMAX 89 103 91).</p>

- <p align="justify">RESCALE_FORCES caps the force on any atom at 0.042 hartree/bohr (about 50 kcal mol<sup>-1</sup> &Aring;<sup>-1</sup>), so that a bad contact in the starting structure does not throw atoms across the box in the first steps.</p>

- <p align="justify">The convergence criteria are loose (forces of 10<sup>-2</sup> hartree/bohr): the goal is to remove the strain, not to reach a minimum. The trajectory is written every 100 steps (EM-pos-1.dcd) and the restart file (EM-1.restart) at the end, which is the starting point of the QM/MM minimization.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Running the minimization </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
mpirun -np 4 cp2k.psmp -i em.inp -o em.out
</pre>

<p align="justify">or with the script that runs all the stages of the tutorial, checking the files that each one needs:</p>

<pre style="color: white; background-color: black;">
./run_cp2k.sh -h
./run_cp2k.sh -s em -n 4                       # this stage only
./run_cp2k.sh -n 64 -e cp2k.psmp               # all stages: em em-me em-ee react smd
./run_cp2k.sh -l srun -n 128 -s "em-ee react"  # inside a Slurm job
./run_cp2k.sh -d                               # print the commands without running them
</pre>

<br/>

- <p align="justify">The cp2k_opt function of cp2k_funcs prints the convergence criteria at every optimization step:</p>

<pre style="color: white; background-color: black;">
source cp2k_funcs
cp2k_opt em.out
</pre>

<br/>

- <p align="justify">Read the beginning of em.out: the number of atoms, kinds and molecules, the cell and the charge of the system. The warning "A non Critical ForceField parameter is missing" is expected (Amber has no Urey-Bradley terms); activate &MM/&PRINT/&FF_INFO to list the missing terms and make sure that no bond, angle or torsion is among them. If CP2K warns that the CRD file has no velocity or box information, the warning is harmless: the cell is given in the input and the minimization needs no velocities.</p>

- <p align="justify">Do not write an inline comment after PARM_FILE_NAME or CONN_FILE_NAME: CP2K takes it as part of the file name and stops with "Force field file missing".</p>

- <p align="justify">Look at the trajectory in VMD and check that the protein, the oligosaccharide and the ions stay where they should:</p>

<pre style="color: white; background-color: black;">
vmd -parm7 amylase.prmtop -dcd EM-pos-1.dcd
</pre>

<br/>
