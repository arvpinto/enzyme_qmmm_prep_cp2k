<br>
<h2> <p align="center"> <b>VIII - Steered QM/MM MD: breaking the glycosidic bond</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
REACT-NVT-EE-GPW-1.restart(.bin) + .wfn  ->  smd-nvt-ee_gpw.inp (moving restraint on the CV, 2 ps)  ->  SMD-NVT-EE-GPW-COLVAR.metadynLog + SMD-NVT-EE-GPW-pos-1.dcd
</pre>

---

<br/>

<p align="justify">In steered MD a harmonic restraint on a collective variable (CV) is moved at constant speed, and the system is dragged from the reactant to the product along it. A single fast pull does not give a free energy, but it shows whether the QM/MM model can describe the reaction, which bonds form and break and in which order, and it gives the structures from which umbrella sampling windows can be started.</p>

<br/>

<h3> <p align="center"><b> The collective variable </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
  &SUBSYS
    ...
    &COLVAR 1               ! CV1
      &DISTANCE_FUNCTION    ! CV = d(1-2) + COEFFICIENT*d(3-4)
        ATOMS 7699 7700 3607 7699 ! CV1 = d(O4-C1, glycosidic bond) - d(HE2(Glu232)-O4, proton transfer)
        COEFFICIENT -1.
      &END DISTANCE_FUNCTION
    &END COLVAR
  &END SUBSYS
</pre>

<br/>

| Serial | Atom | Role |
| :--- | :--- | :--- |
| 7699 | O4 of 4GA500 | glycosidic oxygen: loses C1, gains the proton |
| 7700 | C1 of 4GA501 | anomeric carbon of the -1 sugar |
| 3607 | HE2 of GLH232 | proton of the general acid |

<br/>

<p align="justify">DISTANCE_FUNCTION combines two distances, CV = d(O4-C1) - d(HE2-O4). It is negative in the reactant (bond formed, proton on Glu232: -0.45 &Aring; in the starting structure) and positive when the glycosidic bond is broken and O4 is protonated (e.g. 3.0 - 1.0 = +2.0 &Aring;). A single CV that combines the bond that breaks and the bond that forms lets the system choose the timing of the two events. The attack of the nucleophile, Asp196 OD1 on C1 (3.1 &Aring; in the starting structure), is not part of the CV: follow it along the run to see whether the covalent glycosyl-enzyme intermediate forms or whether the -1 sugar stays as an oxocarbenium-like species.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> The moving restraint </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
&MOTION
  &CONSTRAINT
    &COLLECTIVE
      COLVAR 1
      INTERMOLECULAR TRUE
      TARGET [angstrom] -2.          ! value at t=0; set it close to the CV at the end of the NVT run
      TARGET_GROWTH [angstrom*fs^-1] .002  ! STEPS = (TARGET_LIMIT - TARGET)/TARGET_GROWTH/TIMESTEP
      TARGET_LIMIT [angstrom] 2.
      &RESTRAINT
        K [angstrom^-2*kcalmol] 100.0  ! restraint energy K*(CV-TARGET)^2
      &END RESTRAINT
    &END COLLECTIVE
  &END CONSTRAINT
  &FREE_ENERGY
    METHOD METADYN     ! no hills (DO_HILLS FALSE): only used to print the CV to SMD-NVT-EE-GPW-COLVAR.metadynLog
    &METADYN
      DO_HILLS FALSE
      &METAVAR
        COLVAR 1
        SCALE .5
      &END METAVAR
      &PRINT
        &COLVAR
          COMMON_ITERATION_LEVELS 3
          &EACH
            MD 10
          &END EACH
        &END COLVAR
      &END PRINT
    &END METADYN
  &END FREE_ENERGY
  &MD
    ...
    STEPS 2000
    TIMESTEP 1.
    ...
  &END MD
  ...
&END MOTION

&EXT_RESTART
  BINARY_RESTART_FILE_NAME REACT-NVT-EE-GPW-1.restart.bin
  RESTART_CELL
  RESTART_FILE_NAME REACT-NVT-EE-GPW-1.restart
  RESTART_POS
  RESTART_QMMM
  RESTART_COUNTERS FALSE
  RESTART_THERMOSTAT
&END EXT_RESTART
</pre>

<br/>

- <p align="justify">&COLLECTIVE with a &RESTRAINT section turns the constraint into a harmonic restraint, E = K (CV - TARGET)<sup>2</sup> (CP2K does not use the factor 1/2), with K = 100 kcal mol<sup>-1</sup> &Aring;<sup>-2</sup>. TARGET starts at -2.0 &Aring; and grows by 0.002 &Aring; per fs until TARGET_LIMIT (+2.0 &Aring;), which takes the 2000 steps of 1 fs of the run. INTERMOLECULAR TRUE is needed because the atoms of the CV belong to different molecules (the protein and the oligosaccharide).</p>

- <p align="justify">TARGET should be close to the value of the CV at the start of the run. With TARGET -2.0 and a CV around -0.5 &Aring; (Section VII), the restraint first pushes the CV back towards -2 with a force of 2K x 1.5 = 300 kcal mol<sup>-1</sup> &Aring;<sup>-1</sup>, and the first part of the run is spent undoing that. Set TARGET to the CV at the end of the NVT run and adjust STEPS to (TARGET_LIMIT - TARGET) / TARGET_GROWTH, e.g. 1250 steps from -0.5 to +2.0 &Aring;.</p>

- <p align="justify">The pulling speed (2 &Aring;/ps) and K decide how far from equilibrium the run is. A stiff restraint keeps the CV close to the target; a slower pull leaves more time for the environment to respond, at a higher cost. For the qualitative picture that this run gives, the values of the input are a reasonable start.</p>

- <p align="justify">The &FREE_ENERGY section with DO_HILLS FALSE adds no bias: it is used only to print the CV every 10 MD steps to SMD-NVT-EE-GPW-COLVAR.metadynLog (column 1 is the time in fs, column 2 the CV in bohr).</p>

- <p align="justify">The run continues the NVT simulation: positions, velocities (from the .bin file) and the thermostat are read from the restart files of Section VII, and the wavefunction from REACT-NVT-EE-GPW-RESTART.wfn. RESTART_COUNTERS FALSE starts the step counter and the time at 0, which is what TARGET_GROWTH refers to.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Running the steered MD </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
mpirun -np 128 cp2k.psmp -i smd-nvt-ee_gpw.inp -o smd-nvt-ee_gpw.out
</pre>

<pre style="color: white; background-color: black;">
./run_cp2k.sh -n 128 -s smd
</pre>

<br/>

- <p align="justify">During the run, the CV can be followed with the cp2k_colvar function, which converts it to &Aring;ngstr&ouml;m:</p>

<pre style="color: white; background-color: black;">
source cp2k_funcs
cp2k_colvar SMD-NVT-EE-GPW-COLVAR.metadynLog | tail
</pre>

<br/>

- <p align="justify">As the bonds rearrange, the SCF needs more iterations: check that it keeps converging (grep "SCF run NOT converged" smd-nvt-ee_gpw.out), and stop the run if it does not.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Analysis </b></p> </h3>

<br/>

<p align="justify">plot_smd.py plots the CV against the moving target, the work done by the restraint and, if the .ener file is given, the potential energy, in panels that share the time axis. The plot is shown on screen and saved at 600 dpi:</p>

<pre style="color: white; background-color: black;">
python3 plot_smd.py -c SMD-NVT-EE-GPW-COLVAR.metadynLog -e SMD-NVT-EE-GPW-1.ener \
    --target -2.0 --growth 0.002 --limit 2.0 -k 100 -o smd.png --csv smd.csv
</pre>

<br/>

- <p align="justify">--target, --growth, --limit and -k must match the &COLLECTIVE section of the input. The script also prints the CV at the start and at the end, the largest deviation from the target and the maximum of the work. --csv saves time, CV, target, restraint energy and work to a text file, and --no-show only saves the figure.</p>

- <p align="justify">The work is W(t) = &Sigma; -2K (CV - TARGET) &Delta;TARGET. Its maximum is an upper bound of the free energy barrier along the CV, inflated by the dissipation of a fast, single pull: use it to compare setups, not as a barrier. A potential of mean force needs equilibrium sampling, e.g. umbrella sampling with windows started from frames of this run.</p>

- <p align="justify">The CV going smoothly through zero, with the deviation from the target staying small, means the restraint is in control. A CV that lags behind and then jumps suggests a barrier that the pull crossed abruptly, or a coordinate that is missing from the CV.</p>

<br/>

<p align="justify">The bonds that form and break are best seen in VMD. vmd_cp2k-repr.tcl loads the topology and the trajectories, draws the QM region with dynamic bonds (bonds are drawn only when the atoms are within the bonding distance, so they appear and disappear along the trajectory) and labels the distances given with -bonds, whose values along the trajectory can be plotted from Graphics > Labels > Bonds > Graph:</p>

<pre style="color: white; background-color: black;">
vmd -e vmd_cp2k-repr.tcl -args amylase_ee.prmtop REACT-NVT-EE-GPW-pos-1.dcd SMD-NVT-EE-GPW-pos-1.dcd \
    -qminc forceeval_qmmm.cp2k.inc -bonds "7699:7700 3607:7699 3026:7700"
</pre>

<br/>

| Label | Distance | Expected along the run |
| :--- | :--- | :--- |
| 7699:7700 | O4(500)-C1(501), glycosidic bond | grows from about 1.5 &Aring;: the bond breaks |
| 3607:7699 | HE2(Glu232)-O4(500) | shrinks to about 1.0 &Aring;: O4 is protonated |
| 3026:7700 | OD1(Asp196)-C1(501), nucleophilic attack | shrinks if the glycosyl-enzyme intermediate forms |

<br/>

- <p align="justify">Without -qminc, the script asks for the QM selection in the terminal (the VMD selection of Section II, or the output of cp2k_qmidx). -sel and -selfile give it on the command line, and -step loads every n-th frame of long trajectories.</p>

- <p align="justify">Watch also the rest of the active site: the conformation of the -1 sugar (the ring flattens towards a half-chair as C1 becomes oxocarbenium-like), the hydrogen bonds of Asp299 with the O2 and O3 hydroxyl groups, and the water of the QM region.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Where to go next </b></p> </h3>

<br/>

- <p align="justify">Repeat the pull from different frames of the NVT run, or in the reverse direction, to see whether the same events happen in the same order.</p>

- <p align="justify">For a free energy profile, start umbrella sampling windows along the CV from frames of the steered MD (fixed TARGET, TARGET_GROWTH 0, the same restraint form) and combine them with WHAM or umbrella integration.</p>

- <p align="justify">If the nucleophile does not attack within the pull, add its distance to the CV or use two CVs, and consider the deglycosylation step, in which a water molecule activated by Glu232 hydrolyses the glycosyl-enzyme intermediate.</p>

- <p align="justify">The <a href="https://arvpinto.github.io/enzyme_neb_cp2k/" target="_blank">TS search tutorial</a> shows how to locate and characterize transition states with static QM/MM calculations (CI-NEB and the dimer method) on structures like these.</p>

<br/>
