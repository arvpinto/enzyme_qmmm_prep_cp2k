<br>
<h2> <p align="center"> <b>VIII - Steered QM/MM MD and the Jarzynski work profile</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
QMMM_MD_EQ-1.restart(.bin) + .wfn  ->  qmmm_md_smd.inp (moving restraint on the CV, 2 ps)  ->  QMMM_MD_SMD-COLVAR.metadynLog  ->  smd_jarzynski.py
</pre>

---

<br/>

<p align="justify">In steered MD a harmonic restraint on a collective variable (CV) is moved at constant speed, and the system is dragged from the reactant to the product along it. A single fast pull does not give a free energy, but it shows whether the QM/MM model can describe the reaction, which bonds form and break and in which order, and it gives the structures from which umbrella sampling windows can be started. The run and its analysis are in <a href="https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/tree/main/qmmm_md_smd" target="_blank">qmmm_md_smd</a>.</p>

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
| 7699 | O4 of 4GA500 | glycosidic oxygen: loses C1, gains a proton |
| 7700 | C1 of 4GA501 | anomeric carbon of the -1 sugar |
| 3607 | HE2 of GLH232 | proton of the general acid |

<br/>

<p align="justify">DISTANCE_FUNCTION combines two distances, CV = d(O4-C1) - d(HE2-O4). It is negative in the reactant (bond formed, proton on Glu232: -2.25 &Aring; at the end of Section VII) and positive when the glycosidic bond is broken and O4 is protonated (e.g. 3.0 - 1.0 = +2.0 &Aring;). A single CV that combines the bond that breaks and the bond that forms lets the system choose the timing of the two events. The attack of the nucleophile, Asp196 OD1 on C1 (3.1 &Aring; at the start), is not part of the CV: we follow it along the run.</p>

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
    METHOD METADYN     ! no hills (DO_HILLS FALSE): only used to print the CV to QMMM_MD_SMD-COLVAR.metadynLog
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
  BINARY_RESTART_FILE_NAME QMMM_MD_EQ-1.restart.bin
  RESTART_CELL
  RESTART_FILE_NAME QMMM_MD_EQ-1.restart
  RESTART_POS
  RESTART_QMMM
  RESTART_COUNTERS FALSE
  RESTART_THERMOSTAT
&END EXT_RESTART
</pre>

<br/>

- <p align="justify">&COLLECTIVE with a &RESTRAINT section turns the constraint into a harmonic restraint, U = K (CV - TARGET)<sup>2</sup> (CP2K does not use the factor 1/2), with K = 100 kcal mol<sup>-1</sup> &Aring;<sup>-2</sup>. TARGET starts at -2.0 &Aring; and grows by 0.002 &Aring; per fs (2 &Aring;/ps) until TARGET_LIMIT (+2.0 &Aring;), which takes the 2000 steps of 1 fs of the run. INTERMOLECULAR TRUE is needed because the atoms of the CV belong to different molecules (the protein and the oligosaccharide).</p>

- <p align="justify">TARGET should be close to the value of the CV at the start of the run, otherwise the first part of the pull is spent bringing the CV to the restraint. Here the CV starts at -2.25 &Aring;, 0.25 &Aring; below TARGET, and reaches the target within the first 20 fs. If you change the system or the equilibration, measure the CV at the end of the NVT run (Section VII) and set TARGET and STEPS accordingly.</p>

- <p align="justify">The &FREE_ENERGY section with DO_HILLS FALSE adds no bias: it is used only to print the CV every 10 MD steps to QMMM_MD_SMD-COLVAR.metadynLog (column 1 is the time in fs, column 2 the CV in bohr).</p>

- <p align="justify">The run continues the NVT simulation: positions, velocities (from the .bin file) and the thermostat are read from the restart files of Section VII, and the wavefunction from QMMM_MD_EQ-RESTART.wfn. RESTART_COUNTERS FALSE starts the step counter and the time at 0, which is what TARGET_GROWTH refers to.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Running the steered MD </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
cp -r ~/enzyme_qmmm_prep_cp2k/qmmm_md_smd ~/scratch/ && cd ~/scratch/qmmm_md_smd
cp ~/enzyme_qmmm_prep_cp2k/system_prep/{amylase_ee.prmtop,forceeval_qmmm.cp2k.inc} .
cp ~/scratch/qmmm_md_eq/{QMMM_MD_EQ-1.restart,QMMM_MD_EQ-1.restart.bin,QMMM_MD_EQ-RESTART.wfn} .
mpirun -np 128 cp2k.psmp -i qmmm_md_smd.inp -o qmmm_md_smd.out
</pre>

<pre style="color: white; background-color: black;">
./scripts/run_cp2k.sh -w ~/scratch/amylase -s qmmm_md_smd -n 128
</pre>

<br/>

<p align="justify">The example took 5.8 h on 128 MPI ranks, with every SCF converged. During the run, follow the CV with the cp2k_colvar function, which converts it to &Aring;ngstr&ouml;m:</p>

<pre style="color: white; background-color: black;">
source cp2k_funcs
cp2k_colvar QMMM_MD_SMD-COLVAR.metadynLog | sed -n '1,3p;$p'
#   time[fs] CV[angstrom]
       0.000     -2.25313
      10.000     -2.05654
    2000.000      1.95699
</pre>

<br/>

<p align="justify">After the first 50 fs the CV stays on average 0.07 &Aring; from the moving target: the restraint is in control throughout the pull. The temperature of the system stays at 299.8 &plusmn; 0.8 K, while the QM region again runs hotter (381 K on average, cp2k_eqmmm).</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Watching the bonds form and break </b></p> </h3>

<br/>

<p align="justify">Download QMMM_MD_SMD-pos-1.dcd from the <a href="https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/releases/latest" target="_blank">latest release</a> into qmmm_md_smd. vmd_cp2k-repr.tcl draws the QM region with dynamic bonds, which are drawn only while two atoms are within the bonding distance, so bonds appear and disappear along the trajectory. The labels given with -bonds can be plotted from Graphics > Labels > Bonds > Graph; here they follow the two distances of the CV, the nucleophile and the water of the QM region:</p>

<pre style="color: white; background-color: black;">
cd qmmm_md_smd
wget https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/releases/latest/download/QMMM_MD_SMD-pos-1.dcd
vmd -e ../scripts/vmd_cp2k-repr.tcl -args ../system_prep/amylase_ee.prmtop QMMM_MD_SMD-pos-1.dcd \
    -qminc ../system_prep/forceeval_qmmm.cp2k.inc \
    -bonds "7699:7700 3607:7699 3026:7700 3607:7783 7785:7699"
</pre>

<br/>

| Label | Distance | Start of the pull | End of the pull |
| :--- | :--- | ---: | ---: |
| 7699:7700 | O4(500)-C1(501), glycosidic bond | 1.48 &Aring; | 4.74 &Aring; |
| 3607:7699 | HE2(Glu232)-O4(500) | 3.73 &Aring; | 2.78 &Aring; |
| 3026:7700 | OD1(Asp196)-C1(501), nucleophilic attack | 3.13 &Aring; | 1.53 &Aring; |
| 3607:7783 | HE2(Glu232)-O(WAT511) | 1.68 &Aring; | 0.98 &Aring; |
| 7785:7699 | H2(WAT511)-O4(500) | 2.07 &Aring; | 0.99 &Aring; |

<p align="center"><i>Start: QMMM_MD_EQ-1.restart.bin; end: QMMM_MD_SMD-1.restart.bin.</i></p>

<br/>

<p align="justify">At the end of the pull the glycosylation step is complete, but not exactly as the CV describes it:</p>

- <p align="justify">the glycosidic bond is broken, and Asp196 OD1 is bonded to C1 (1.53 &Aring;): the covalent glycosyl-enzyme intermediate has formed, although the nucleophile is not in the CV. C1-O5 of the -1 sugar is 1.34 &Aring;, shorter than at the start (1.39 &Aring;);</p>

- <p align="justify">O4 is protonated, but not by HE2 of Glu232. At the start of the pull the water of the QM region already bridges the two (Glu232 HE2...O<sub>W</sub> 1.68 &Aring;, H2<sub>W</sub>...O4 2.07 &Aring;); at the end HE2 belongs to the water (0.98 &Aring;, Glu232 OE2-HE2 1.85 &Aring;) and one of the original water protons belongs to O4 (0.99 &Aring;). The proton went from Glu232 to O4 through WAT511, a Grotthuss-like relay.</p>

<p align="justify">Play the trajectory to see when each event happens (the bond lengths only tell us where the run ended). The relay has a consequence for the CV: it follows HE2, which never reaches O4 (2.78 &Aring; at the end), so once the proton has moved through the water the second term of the CV no longer describes the proton transfer, and the CV grows by stretching O4 away from C1 (4.74 &Aring;). This is a common pitfall of CVs defined on specific hydrogen atoms; a CV that does not depend on which proton arrives (e.g. a coordination number of O4 by hydrogens) avoids it. Watch also the conformation of the -1 sugar ring and the hydrogen bonds of Asp299 with its O2 and O3 hydroxyls.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Following the pull in time </b></p> </h3>

<br/>

<p align="justify">plot_smd.py puts on a common time axis what you see in VMD: the CV and the restraint centre, the work done by the restraint and, if the trajectory is given, the distances of your choice. Read the panels next to the animation (the DCD has one frame every 5 fs, so t = 0.72 ps, the first maximum of the work, is around frame 144) to see what the system is doing when the work rises or drops. The restraint parameters, the time step and the trajectory stride are read from the input echo in qmmm_md_smd.out:</p>

<pre style="color: white; background-color: black;">
python3 ../scripts/plot_smd.py -h
python3 ../scripts/plot_smd.py -l QMMM_MD_SMD-COLVAR.metadynLog -o qmmm_md_smd.out
</pre>

<pre style="color: white; background-color: black;">
Restraint: TARGET -2.000 Å, TARGET_LIMIT 2.000 Å, K 100.00 kcal/mol/Å^2 (U = K (s-lambda)^2), speed 2.000 Å/ps; time step 1 fs
QMMM_MD_SMD-COLVAR.metadynLog: 201 points, t = 0.0 to 2000.0 fs
  CV     start   -2.253  end    1.957  min   -2.253  max    1.974 Å
  |CV - lambda|  mean 0.066 Å, max 0.253 Å at t = 0.0 fs
  work   end    26.09 kcal/mol
  work   local max    15.19 kcal/mol at t =   720.0 fs, lambda =  -0.56 Å
  work   local min    11.71 kcal/mol at t =   990.0 fs, lambda =  -0.02 Å
  work   local max    12.72 kcal/mol at t =  1040.0 fs, lambda =   0.08 Å
  work   local min    10.94 kcal/mol at t =  1210.0 fs, lambda =   0.42 Å
wrote plot_smd.dat
wrote plot_smd.png (600 dpi)
</pre>

<br/>

<div align="center">
    <img src="assets/img/plot_smd.png" width="80%">
</div>

<br/>

- <p align="justify">The work is integrated as in smd_jarzynski.py (trapezoidal rule on the COLVAR points), so both scripts give the same numbers. The local maxima and minima are the points that are the highest or lowest within &plusmn;150 fs (--window).</p>

- <p align="justify">With the trajectory, the distances of the table above are added as a third panel. The pairs use the serial numbers of the -bonds option of vmd_cp2k-repr.tcl, and their values in every frame are written to plot_smd_dist.dat:</p>

<pre style="color: white; background-color: black;">
python3 ../scripts/plot_smd.py -l QMMM_MD_SMD-COLVAR.metadynLog -o qmmm_md_smd.out \
    --dcd QMMM_MD_SMD-pos-1.dcd --pairs "7699:7700 3607:7699 3026:7700 3607:7783 7785:7699" \
    --pair-labels "O4-C1,HE2-O4,OD1-C1,HE2-Ow,H2w-O4"
</pre>

<p align="justify">This panel shows the order of the events directly: when the glycosidic bond breaks, when OD1 reaches C1, and when each proton of the relay moves.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Work and free energy: the Jarzynski equality </b></p> </h3>

<br/>

<p align="justify">The work done by the moving restraint is W(t) = &int; &part;U/&part;&lambda; d&lambda;, with &lambda;(t) = TARGET + TARGET_GROWTH t and &part;U/&part;&lambda; = -2K (CV - &lambda;). The Jarzynski equality relates the work of N independent pulls to the free energy along &lambda;, F(&lambda;) = -kT ln &lt;exp(-W(&lambda;)/kT)&gt;. smd_jarzynski.py integrates the work from the metadynLog files, reads TARGET, TARGET_GROWTH, TARGET_LIMIT, K and the temperature from the echo of the input in the CP2K output (or from the command line), and computes the Jarzynski average, the second-order cumulant expansion, F = &lt;W&gt; - var(W)/2kT, and their bootstrap errors when there are several replicas:</p>

<pre style="color: white; background-color: black;">
python3 ../scripts/smd_jarzynski.py -h
python3 ../scripts/smd_jarzynski.py -l QMMM_MD_SMD-COLVAR.metadynLog -o qmmm_md_smd.out
</pre>

<pre style="color: white; background-color: black;">
# TARGET = -3.779452  TARGET_GROWTH = 9.142058e-05  TARGET_LIMIT = 3.779452265771286  K = 4.462539e-02  T = 300.0 K
# K = 100.00 kcal/mol/Å^2 as used by CP2K (U = K (s-lambda)^2); equivalent k = 200.00 in the U = k/2 (s-lambda)^2 convention
# pulling speed = 2.000 Å/ps
# QMMM_MD_SMD-COLVAR.metadynLog: 201 points, final W = 26.09 kcal/mol
# NOTE: single trajectory -> F(lambda) = W(lambda), an upper bound that includes dissipation
# wrote smd_jarzynski.dat
# wrote smd_jarzynski.png (600 dpi)
</pre>

<br/>

<div align="center">
    <img src="assets/img/smd_jarzynski.png" width="100%">
</div>

<br/>

- <p align="justify">The parameters are read in CP2K internal units (bohr and hartree); the script prints them back in &Aring; and kcal/mol as a check. The profile is plotted against &lambda;, the restraint center, which is the proper variable of the Jarzynski equality (--xaxis cv plots it against the CV, averaged over the replicas, instead). smd_jarzynski.dat has the time, &lambda;, the CV, the work, both estimators and their errors at every point of the COLVAR file; the figure is shown on screen and saved at 600 dpi (--no-show only saves it, --prefix changes the file names).</p>

- <p align="justify">With a single trajectory both estimators reduce to the work itself, an upper bound of the free energy that includes the dissipation of a fast pull. In the example, the work rises to 15.2 kcal/mol at &lambda; = -0.56 &Aring; (t = 0.72 ps), drops to 10.9 kcal/mol at &lambda; = +0.42 &Aring; (t = 1.21 ps) and then rises again to 26.1 kcal/mol at the end, as the restraint keeps stretching O4 away from C1 after the relay (see above). The first maximum should not be read as the barrier of glycosylation: the work of a single fast pull includes dissipation, so it overestimates the free energy at every &lambda;.</p>

- <p align="justify">For an actual free energy estimate, run several pulls with the same protocol from different starting points (e.g. frames of a longer NVT run, or new velocities), give all their metadynLog files to the script, and check that the Jarzynski and cumulant estimates agree within their bootstrap errors; slower pulls need fewer replicas. Umbrella sampling with windows started from frames of the steered MD is the usual alternative.</p>

<pre style="color: white; background-color: black;">
python3 ../scripts/smd_jarzynski.py -l rep*/QMMM_MD_SMD-COLVAR.metadynLog -o rep1/qmmm_md_smd.out --nboot 500 --prefix smd_replicas
</pre>

<br/>

---

<br/>

<h3> <p align="center"><b> Where to go next </b></p> </h3>

<br/>

- <p align="justify">Repeat the pull from different frames of the NVT run, or more slowly, and combine the replicas with smd_jarzynski.py.</p>

- <p align="justify">Redefine the CV so that it does not depend on which proton reaches O4 (e.g. with a coordination number), or add the nucleophile distance as a second CV, and compare which events happen and in which order.</p>

- <p align="justify">Start umbrella sampling windows along the CV from frames of the steered MD (fixed TARGET, TARGET_GROWTH 0, the same restraint form) and combine them with WHAM or umbrella integration.</p>

- <p align="justify">Study the deglycosylation step, in which a water molecule activated by Glu232 hydrolyses the glycosyl-enzyme intermediate formed here.</p>

- <p align="justify">The <a href="https://arvpinto.github.io/enzyme_neb_cp2k/" target="_blank">TS search tutorial</a> shows how to locate and characterize transition states with static QM/MM calculations (CI-NEB and the dimer method) on structures like these.</p>

<br/>
