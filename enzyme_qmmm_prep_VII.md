<br>
<h2> <p align="center"> <b>VII - QM/MM MD of the reactant</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
EM-EE-GPW-1.restart + EM-EE-GPW-RESTART.wfn  ->  react_nvt-ee_gpw.inp (NVT, 300 K, 2 ps)  ->  REACT-NVT-EE-GPW-1.restart(.bin) + REACT-NVT-EE-GPW-RESTART.wfn
</pre>

---

<br/>

<p align="justify">Before pulling the system along the reaction coordinate, the reactant is equilibrated with QM/MM MD at 300 K, so that the steered MD starts from a thermalized structure with velocities, and so that we know where the collective variable sits in the reactant state.</p>

<br/>

<h3> <p align="center"><b> The input </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
&GLOBAL
  ECHO_INPUT
  PRINT_LEVEL LOW
  PROJECT REACT-NVT-EE-GPW
  RUN_TYPE MD
  SAVE_MEM
&END GLOBAL
...
&MOTION
  &MD
    DISPLACEMENT_TOL 5.29177209E+001
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
  RESTART_FILE_NAME EM-EE-GPW-1.restart
  RESTART_POS
  RESTART_QMMM
  RESTART_COUNTERS FALSE
  RESTART_THERMOSTAT FALSE
&END EXT_RESTART
</pre>

<br/>

- <p align="justify">The &FORCE_EVAL section is that of em-ee_gpw.inp (electrostatic embedding with amylase_ee.prmtop), with the wavefunction of the minimization as the SCF guess. In MD the wavefunction of the previous step is a good guess, so the SCF loops are shorter (MAX_SCF 30 and 10 outer loops).</p>

- <p align="justify">2000 steps of 1 fs at 300 K in the NVT ensemble, with the canonical sampling through velocity rescaling (CSVR) thermostat and a time constant of 100 fs. With REGION DEFINED, each &DEFINE_REGION section is one thermostat region; the single section of the input holds the QM atoms (QM_SUBSYS ATOMIC) and the MM atoms (MM_SUBSYS ATOMIC), so one thermostat acts on the whole system. To thermostat the QM region on its own (its temperature fluctuates much more than that of the 68 000 MM atoms), split it into two &DEFINE_REGION sections, one with QM_SUBSYS ATOMIC and one with MM_SUBSYS ATOMIC. The minimization restart has no velocities, so the run starts from velocities generated at TEMPERATURE.</p>

- <p align="justify">The restart file is written every 100 steps (SPLIT_RESTART_FILE writes positions and velocities to the binary REACT-NVT-EE-GPW-1.restart.bin, which keeps the text restart file small), the trajectory every 5 steps (REACT-NVT-EE-GPW-pos-1.dcd) and the energies of every step to REACT-NVT-EE-GPW-1.ener.</p>

- <p align="justify">2 ps is enough to thermalize the region around the QM atoms after the minimization and to start the steered MD from a dynamic structure, not to sample the reactant state. For a free energy study, run longer and/or several replicas from different initial velocities.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Running and checking the MD </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
mpirun -np 128 cp2k.psmp -i react_nvt-ee_gpw.inp -o react_nvt-ee_gpw.out
</pre>

<pre style="color: white; background-color: black;">
./run_cp2k.sh -n 128 -s react
</pre>

<br/>

- <p align="justify">The cp2k_ener function collects the MD information of every step (time, energies, temperatures and CPU time) in a table, and cp2k_eqmmm the QM energy and the temperature of the QM region:</p>

<pre style="color: white; background-color: black;">
source cp2k_funcs
cp2k_ener react_nvt-ee_gpw.out      # -> react_nvt-ee_gpw.ener.csv
cp2k_eqmmm react_nvt-ee_gpw.out     # -> react_nvt-ee_gpw.eqmmm.csv
</pre>

<br/>

- <p align="justify">Check that the temperature of the system and of the QM region (printed separately in the output) fluctuates around 300 K after the first few hundred femtoseconds, that the conserved quantity has no large drift, that every SCF converges (grep "SCF run NOT converged" react_nvt-ee_gpw.out should print nothing), and that no QM atom gets close to the border of the QM cell (SKIN warnings). The water of the QM region (WAT511) and the sugars must stay in place: look at the trajectory with vmd_cp2k-repr.tcl (Section VIII).</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Where does the reactant sit on the CV? </b></p> </h3>

<br/>

<p align="justify">The CV of the next section is the difference between the glycosidic bond C1(501)-O4(500) and the distance between the acidic proton of Glu232 and the same O4 (serials 7700, 7699 and 3607). Before setting the starting value of the restraint, measure it along the NVT trajectory, e.g. in the VMD console (index = serial - 1):</p>

<pre style="color: white; background-color: black;">
vmd -parm7 amylase_ee.prmtop -dcd REACT-NVT-EE-GPW-pos-1.dcd
</pre>

<pre style="color: white; background-color: black;">
set n [molinfo top get numframes]
for {set f 0} {$f &lt; $n} {incr f 20} {
    set d1 [measure bond {7698 7699} frame $f]
    set d2 [measure bond {3606 7698} frame $f]
    puts [format "frame %4d  d(O4-C1) %5.2f  d(HE2-O4) %5.2f  CV %6.2f" $f $d1 $d2 [expr {$d1-$d2}]]
}
</pre>

<br/>

<p align="justify">In the starting structure (frame_xray.pdb) the glycosidic bond is 1.52 &Aring; and the proton of Glu232 is hydrogen bonded to O4 at 1.97 &Aring;, so the CV is about -0.45 &Aring;. Its value at the end of the NVT run is the natural starting point of the restraint (Section VIII).</p>

<br/>
