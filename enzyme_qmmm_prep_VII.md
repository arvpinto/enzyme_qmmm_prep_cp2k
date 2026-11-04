<br>
<h2> <p align="center"> <b>VII - QM/MM MD of the reactant</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
QMMM_EM_EE-1.restart + QMMM_EM_EE-RESTART.wfn  ->  qmmm_md_eq.inp (NVT, 300 K, 2 ps)  ->  QMMM_MD_EQ-1.restart(.bin) + QMMM_MD_EQ-RESTART.wfn
</pre>

---

<br/>

<p align="justify">Before pulling the system along the reaction coordinate, the reactant is equilibrated with QM/MM MD at 300 K, so that the steered MD starts from a thermalized structure with velocities, and so that we know where the collective variable sits in the reactant state. The run is in <a href="https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/tree/main/qmmm_md_eq" target="_blank">qmmm_md_eq</a>.</p>

<br/>

<h3> <p align="center"><b> The input </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
&GLOBAL
  ECHO_INPUT
  PRINT_LEVEL LOW
  PROJECT QMMM_MD_EQ
  RUN_TYPE MD
  SAVE_MEM
&END GLOBAL
...
&MOTION
  &MD
    ENSEMBLE NVT
    MAX_STEPS 10000
    STEPS 2000
    STEP_START_VAL 0
    TEMPERATURE 300.0
    TEMP_KIND FALSE
    TIMESTEP 1.
    TIME_START_VAL 0.0
    &THERMOSTAT
      REGION DEFINED
      TYPE CSVR
      &CSVR
        TIMECON 100
      &END CSVR
      &DEFINE_REGION
        MM_SUBSYS ATOMIC
        QM_SUBSYS ATOMIC
      &END DEFINE_REGION
    &END THERMOSTAT
  &END MD
  &PRINT
    &RESTART
      &EACH
        MD 100
      &END EACH
      ADD_LAST NUMERIC
      SPLIT_RESTART_FILE TRUE
    &END RESTART
    ...
    &TRAJECTORY
      &EACH
        MD 5
      &END EACH
      ADD_LAST NUMERIC
      FORMAT DCD
    &END TRAJECTORY
  &END PRINT
&END MOTION

&EXT_RESTART
  RESTART_CELL
  RESTART_FILE_NAME QMMM_EM_EE-1.restart
  RESTART_POS
  RESTART_QMMM
  RESTART_COUNTERS FALSE
  RESTART_THERMOSTAT FALSE
&END EXT_RESTART
</pre>

<br/>

- <p align="justify">The &FORCE_EVAL section is that of qmmm_em_ee.inp (electrostatic embedding with amylase_ee.prmtop), with the wavefunction of the minimization as the SCF guess and without &RESCALE_FORCES. In MD the wavefunction of the previous step is a good guess, so the SCF loops are shorter (MAX_SCF 30 and 10 outer loops).</p>

- <p align="justify">2000 steps of 1 fs at 300 K in the NVT ensemble, with the canonical sampling through velocity rescaling (CSVR) thermostat and a time constant of 100 fs. With REGION DEFINED, each &DEFINE_REGION section is one thermostat region; the single section of the input holds the QM atoms (QM_SUBSYS ATOMIC) and the MM atoms (MM_SUBSYS ATOMIC), so one thermostat acts on the whole system. The minimization restart has no velocities, so the run starts from velocities generated at TEMPERATURE.</p>

- <p align="justify">The restart file is written every 100 steps (SPLIT_RESTART_FILE writes positions and velocities to the binary QMMM_MD_EQ-1.restart.bin, which keeps the text restart file small), the trajectory every 5 steps (QMMM_MD_EQ-pos-1.dcd) and the energies of every step to QMMM_MD_EQ-1.ener.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Running and checking the MD </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
cp -r ~/enzyme_qmmm_prep_cp2k/qmmm_md_eq ~/scratch/ && cd ~/scratch/qmmm_md_eq
cp ~/enzyme_qmmm_prep_cp2k/system_prep/{amylase_ee.prmtop,forceeval_qmmm.cp2k.inc} .
cp ~/scratch/qmmm_em/{QMMM_EM_EE-1.restart,QMMM_EM_EE-RESTART.wfn} .
mpirun -np 128 cp2k.psmp -i qmmm_md_eq.inp -o qmmm_md_eq.out
</pre>

<pre style="color: white; background-color: black;">
./scripts/run_cp2k.sh -w ~/scratch/amylase -s qmmm_md_eq -n 128
</pre>

<br/>

<p align="justify">In the example, each MD step took about 10.5 s on 128 MPI ranks (5.9 h for the 2 ps), and every SCF converged (grep "SCF run NOT converged" qmmm_md_eq.out prints nothing). The cp2k_ener function collects the MD information of every step (time, energies, temperatures and CPU time) in a table, and cp2k_eqmmm the QM energy and the temperature of the QM region:</p>

<pre style="color: white; background-color: black;">
source cp2k_funcs
cp2k_ener qmmm_md_eq.out      # -> qmmm_md_eq.ener.csv
cp2k_eqmmm qmmm_md_eq.out     # -> qmmm_md_eq.eqmmm.csv
tail -1 qmmm_md_eq.eqmmm.csv
2000.000000 -448.51471611161202 0.119478476414E+00 0.114309351340E+00 359.316894 343.771383
</pre>

<br/>

- <p align="justify">The temperature of the system fluctuates around 300 K after the first 0.5 ps (296.8 &plusmn; 4.7 K, from QMMM_MD_EQ-1.ener). The conserved quantity changes by about 4200 kcal/mol in the first 0.5 ps, while the minimized structure heats up; after that its net drift is below 100 kcal/mol, with fluctuations within about 1100 kcal/mol (0.2 % of the total).</p>

- <p align="justify">The QM region runs hotter than the rest of the system: its average temperature over the run is 344 K (last column of qmmm_md_eq.eqmmm.csv), and 369 K over the last 1.5 ps. With a single thermostat for 68 000 atoms, the heat released or exchanged in a region of 65 atoms is not removed where it appears. To thermostat the QM region on its own, split the region into two &DEFINE_REGION sections, one with QM_SUBSYS ATOMIC and one with MM_SUBSYS ATOMIC, and compare the QM temperature again.</p>

- <p align="justify">2 ps thermalizes the active site after the minimization and gives the steered MD a dynamic starting point; it does not sample the reactant state. For a free energy study, run longer and/or several replicas from different initial velocities.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Looking at the trajectory, and where the reactant sits on the CV </b></p> </h3>

<br/>

<div class="qm-movie">
  <video src="assets/movies/qmmm_md_eq_qm.mp4" autoplay loop muted playsinline controls preload="metadata" aria-label="QM region during the QM/MM MD at 300 K"></video>
</div>

<p align="center"><i>The QM region along the 2 ps of QM/MM MD (QMMM_MD_EQ-pos-1.dcd, one frame every 5 fs), smoothed with a moving average over 5 frames (25 fs), with the same view and dashes as in Sections V and VI.</i></p>

<br/>

<p align="justify">Download QMMM_MD_EQ-pos-1.dcd from the <a href="https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/releases/latest" target="_blank">latest release</a> into qmmm_md_eq and play it with the QM region drawn with dynamic bonds:</p>

<pre style="color: white; background-color: black;">
cd qmmm_md_eq
wget https://github.com/arvpinto/enzyme_qmmm_prep_cp2k/releases/latest/download/QMMM_MD_EQ-pos-1.dcd
vmd -e ../scripts/vmd_cp2k-repr.tcl -args ../system_prep/amylase_ee.prmtop QMMM_MD_EQ-pos-1.dcd \
    -qminc ../system_prep/forceeval_qmmm.cp2k.inc -bonds "7699:7700 3607:7699 3026:7700"
</pre>

<p align="justify">The CV of the next section is the difference between the glycosidic bond C1(501)-O4(500) and the distance between the acidic proton of Glu232 and the same O4 (serials 7700, 7699 and 3607). Measure it along the trajectory in the VMD console (index = serial - 1), or plot the labels from Graphics > Labels > Bonds > Graph:</p>

<pre style="color: white; background-color: black;">
set n [molinfo top get numframes]
for {set f 0} {$f &lt; $n} {incr f 20} {
    set d1 [measure bond {7698 7699} frame $f]
    set d2 [measure bond {3606 7698} frame $f]
    puts [format "frame %4d  d(O4-C1) %5.2f  d(HE2-O4) %5.2f  CV %6.2f" $f $d1 $d2 [expr {$d1-$d2}]]
}
</pre>

<br/>

<p align="justify">At the end of the run (QMMM_MD_EQ-1.restart.bin), d(O4-C1) = 1.48 &Aring; and d(HE2-O4) = 3.73 &Aring;, so the CV is -2.25 &Aring;: the proton of Glu232 moved further away from the glycosidic oxygen than in the minimized structure (-1.76 &Aring;, Section VI). This value is the first point of the COLVAR file of the steered MD, and it sets where the pull of Section VIII starts.</p>

<br/>
