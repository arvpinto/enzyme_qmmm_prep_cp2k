<br>
<h2> <p align="center"> <b>VI - QM/MM minimization: mechanical, then electrostatic embedding</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
EM-1.restart  ->  em-me_gpw.inp (ECOUPL NONE, amylase.prmtop)  ->  em-ee_gpw.inp (ECOUPL GAUSS, amylase_ee.prmtop)  ->  EM-EE-GPW-1.restart + EM-EE-GPW-RESTART.wfn
</pre>

---

<br/>

<p align="justify">The MM-minimized structure is relaxed with the QM/MM Hamiltonian in two steps. The first uses mechanical embedding, where the SCF starts from scratch and is cheap and robust; its wavefunction is then the starting guess of the second, with electrostatic embedding, which gives the structure and the wavefunction from which the QM/MM MD starts.</p>

<br/>

<h3> <p align="center"><b> The QM part </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
&FORCE_EVAL
  METHOD QMMM
  STRESS_TENSOR ANALYTICAL
  &DFT
    BASIS_SET_FILE_NAME BASIS_SET
    CHARGE -2
    MULTIPLICITY 1
    POTENTIAL_FILE_NAME POTENTIAL
    UKS FALSE
    &MGRID
      COMMENSURATE  ! mandatory to use GEEP
      CUTOFF 300    ! Default value: 2.80000000E+002
    &END MGRID
    &QS
      METHOD GPW
    &END QS
    &SCF
      EPS_SCF 1E-6
      MAX_SCF 60
      SCF_GUESS ATOMIC
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
        SAFE_DIIS TRUE
      &END OT
      &OUTER_SCF
        EPS_SCF 1E-6
        MAX_SCF 30
      &END OUTER_SCF
    &END SCF
    &XC
      DENSITY_CUTOFF 1.0E-10
      GRADIENT_CUTOFF 1.0E-10
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  ...
  &SUBSYS
    ...
    &KIND C
      BASIS_SET DZVP-GTH-PBE
      POTENTIAL GTH-PBE-q4
    &END KIND
    &KIND O
      BASIS_SET DZVP-GTH-PBE
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET DZV-GTH-PBE
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
</pre>

<br/>

- <p align="justify">The QM region is described with PBE, double-zeta Gaussian basis sets (DZVP for C and O, DZV for H) and GTH pseudopotentials, with the GPW method and a plane-wave cutoff of 300 Ry for the density. CHARGE -2 and MULTIPLICITY 1 are those of the QM region (Section II). BASIS_SET and POTENTIAL are the files distributed in the data directory of CP2K.</p>

- <p align="justify">The SCF uses the orbital transformation method (OT) with DIIS and the FULL_SINGLE_INVERSE preconditioner, with an outer loop for robustness. In the mechanical embedding run the guess is ATOMIC; every later run reads the wavefunction of the previous one (SCF_GUESS RESTART and WFN_RESTART_FILE_NAME), which saves many SCF steps.</p>

- <p align="justify">COMMENSURATE makes the multigrids of the QM cell commensurate with each other, which the GEEP treatment of the QM/MM electrostatics needs.</p>

- <p align="justify">PBE has no dispersion, so it is missing within the QM region (the QM-MM van der Waals interactions come from the Lennard-Jones terms of the force field). For production work consider PBE-D3 (&XC/&VDW_POTENTIAL) or a hybrid functional, and test the basis set and the cutoff on the energies you are interested in.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> The &QMMM section </b></p> </h3>

<br/>

<p align="justify">Mechanical embedding (em-me_gpw.inp, with amylase.prmtop):</p>

<pre style="color: white; background-color: black;">
  &QMMM
    CENTER EVERY_STEP
    ECOUPL NONE
    USE_GEEP_LIB 6
    @INCLUDE forceeval_qmmm.cp2k.inc
  &END QMMM
</pre>

<p align="justify">Electrostatic embedding (em-ee_gpw.inp, react_nvt-ee_gpw.inp and smd-nvt-ee_gpw.inp, with amylase_ee.prmtop):</p>

<pre style="color: white; background-color: black;">
  &QMMM
    CENTER EVERY_STEP
    EPS_MM_RSPACE 1.0E-10
    ECOUPL GAUSS
    MM_POTENTIAL_FILE_NAME MM_POTENTIAL
    PARALLEL_SCHEME ATOM
    SPHERICAL_CUTOFF -5.29177209E-001 0.00000000E+000
    USE_GEEP_LIB 6
    @INCLUDE forceeval_qmmm.cp2k.inc
  &END QMMM
</pre>

<br/>

- <p align="justify">ECOUPL GAUSS describes each MM point charge as a Gaussian, so that the potential it creates on the QM grid has no singularity, and evaluates it with the Gaussian Expansion of the Electrostatic Potential (GEEP): USE_GEEP_LIB 6 takes an expansion with 6 Gaussians from the library built into CP2K (MM_POTENTIAL, in the data directory, is the file of expansions read when the library is not used). EPS_MM_RSPACE is the accuracy of the real-space part of the expansion and PARALLEL_SCHEME ATOM distributes the MM atoms over the MPI tasks.</p>

- <p align="justify">A negative SPHERICAL_CUTOFF switches the cutoff off: every MM charge acts on the QM density. The value written here (-0.529 &Aring;) is -1 bohr, as CP2K prints it in a restart file.</p>

- <p align="justify">The rest of each input (the &MM section, &RESCALE_FORCES and the &TOPOLOGY of &SUBSYS) is that of em.inp, without COORD_FILE_NAME and without the &CELL include: the positions and the cell come from the restart file of the previous step (&EXT_RESTART with RESTART_POS and RESTART_CELL). RESTART_COUNTERS FALSE starts each run from step 0.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Running and checking the minimizations </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
mpirun -np 64 cp2k.psmp -i em-me_gpw.inp -o em-me_gpw.out
mpirun -np 64 cp2k.psmp -i em-ee_gpw.inp -o em-ee_gpw.out
</pre>

<pre style="color: white; background-color: black;">
./run_cp2k.sh -n 64 -s "em-me em-ee"
</pre>

<br/>

<p align="justify">Before the optimization goes far, check the first SCF of each run, in this order:</p>

<pre style="color: white; background-color: black;">
grep -E "Total charge of the classical|Number of electrons|SCF run converged|outer SCF loop|SKIN|Total FORCE_EVAL" em-ee_gpw.out | head
</pre>

- <p align="justify">the charge of the classical system: +2 with the corrected topology (electrostatic embedding, Section IV);</p>
- <p align="justify">the number of electrons: 216 for this QM region (18 C x 4 + 18 O x 6 + 29 H + 5 capping H + 2 for the charge of -2);</p>
- <p align="justify">the convergence of the SCF, and no warnings about QM atoms within the SKIN of the QM cell;</p>
- <p align="justify">the energy and the convergence criteria along the optimization (cp2k_opt em-me_gpw.out).</p>

<br/>

- <p align="justify">Both inputs stop after 100 optimization steps (MAX_ITER). The mechanical embedding run has looser criteria (MAX_FORCE 1.2 x 10<sup>-3</sup> hartree/bohr) than the electrostatic embedding run (the CP2K defaults, 4.5 x 10<sup>-4</sup>). If a run reaches MAX_ITER without converging, continue it from its own restart file (change RESTART_FILE_NAME and WFN_RESTART_FILE_NAME of a copy of the input) before going on.</p>

- <p align="justify">The energies of the two embedding schemes cannot be compared with each other. What matters is that each run converges and that the structure stays chemically sensible: look at the trajectories (written at every optimization step) with the representation script of Section VIII, and check the capping O-H and C-H distances and the hydrogen bonds of the active site.</p>

<pre style="color: white; background-color: black;">
vmd -e vmd_cp2k-repr.tcl -args amylase_ee.prmtop EM-ME-GPW-pos-1.dcd EM-EE-GPW-pos-1.dcd -qminc forceeval_qmmm.cp2k.inc
</pre>

<br/>
