#!/usr/bin/env vmd
# vmd_cp2k-qmmm.tcl
#
# Writes the CP2K include files of a QM/MM model from a VMD selection of the QM layer:
#   forceeval_subsys.cp2k.inc  &CELL of the MM system (FORCE_EVAL > SUBSYS)
#   forceeval_qmmm.cp2k.inc    &CELL of the QM box, &QM_KIND and &LINK sections (FORCE_EVAL > QMMM)
#   qmmm-ee.parmed.in          ParmEd script with the charge correction for electrostatic embedding
# and adds representations of the QM layer and of the QM/MM bonds to the top molecule.
#
# Usage (topology and coordinates with box information, e.g. the tleap output):
#   vmd -parm7 amylase.prmtop -rst7 frame_xray.inpcrd -e vmd_cp2k-qmmm.tcl
#       -> the script asks for the QM selection (VMD syntax) in the terminal
#   vmd -parm7 amylase.prmtop -rst7 frame_xray.inpcrd -e vmd_cp2k-qmmm.tcl -args -selfile qm_sel.vmd
#   vmd -parm7 amylase.prmtop -rst7 frame_xray.inpcrd -e vmd_cp2k-qmmm.tcl -args -sel "resid 196 and ..."
#   (add -dispdev text to run without graphics)
#
# Options (after -args):
#   -sel "<selection>"   QM layer in VMD selection syntax
#   -selfile <file>      file with the QM layer selection (VMD syntax, may span several lines)
#   -pad <A>             padding added to each side of the QM layer to build the QM cell (default 5.0 A)
#   -scale <f>           QMMM_SCALE_FACTOR of the &LINK sections (default 0.5)
#   -lj "<types>"        hydrogen types that get LJ parameters 0.3019 A / 0.047 kcal/mol in the parmed
#                        script (default "HO,ho HW,hw"; use "none" to skip)
#   -overwrite           replace existing output files (by default the script stops if they exist)
#   -quit                exit VMD at the end
#   -h                   print this help

proc qmmm_usage {} {
    puts "Usage: vmd -parm7 <prmtop> -rst7 <rst7/inpcrd> -e vmd_cp2k-qmmm.tcl \[-args \[-sel \"<sel>\" | -selfile <file>\] \[-pad 5.0\] \[-scale 0.5\] \[-lj \"HO,ho HW,hw\"\] \[-overwrite\] \[-quit\]\]"
}

#--> parse options
set qm_layer ""
set qm_pad 5.0
set qmmm_sc 0.5
set lj_types [list "HO,ho" "HW,hw"]
set overwrite 0
set quit_at_end 0
for {set i 0} {$i < [llength $argv]} {incr i} {
    set opt [lindex $argv $i]
    switch -- $opt {
        -sel       { incr i; set qm_layer [lindex $argv $i] }
        -selfile   { incr i; set fh [open [lindex $argv $i] r]; set qm_layer [string trim [join [split [read $fh] "\n"] " "]]; close $fh }
        -pad       { incr i; set qm_pad [lindex $argv $i] }
        -scale     { incr i; set qmmm_sc [lindex $argv $i] }
        -lj        { incr i; set lj_arg [lindex $argv $i]; if {$lj_arg eq "none"} { set lj_types [list] } else { set lj_types $lj_arg } }
        -overwrite { set overwrite 1 }
        -quit      { set quit_at_end 1 }
        -h         { qmmm_usage; return }
        default    { puts "Unknown option $opt"; qmmm_usage; return }
    }
}

if { [molinfo num] == 0 } {
    puts "ERROR: no molecule loaded. Load the topology and coordinates first."
    qmmm_usage
    return
}
display rendermode GLSL

puts "1) MAKE SURE YOUR SYSTEM IS THE TOP MOLECULE!"
set topid [molinfo top get id]
puts "CONTROL STRING: TOP MOLECULE IS $topid!"

if { $qm_layer eq "" } {
    puts "\n2) INTRODUCE SELECTION FOR QM LAYER IN VMD SELECTION FORMAT:"
    set qm_layer [gets stdin]
} else {
    puts "\n2) QM LAYER SELECTION:\n$qm_layer"
}

#--> output files (opened with EXCL unless -overwrite, so that previous files are not lost)
set out_mode [list WRONLY CREAT EXCL]
if { $overwrite } { set out_mode [list WRONLY CREAT TRUNC] }
foreach f {forceeval_subsys.cp2k.inc forceeval_qmmm.cp2k.inc qmmm-ee.parmed.in} {
    if { !$overwrite && [file exists $f] } {
        puts "ERROR: $f already exists. Remove it or use -args -overwrite."
        return
    }
}

#--> set INCLUDE file for CP2K SUBSYS section
set subsys_list [list]

#--> set SYSTEM cell
lappend subsys_list "! CP2K > FORCE_EVAL > SUBSYS > CELL"
##if cell info is not on coordinate file estimate based on min and max system vectors
#set cell_lim [measure minmax [atomselect $topid "noh all"]]
#set cell_size [vecsub [lindex $cell_lim 1] [lindex $cell_lim 0]]
#pbc set $cell_size -molid $topid
#get cell from coordinate file (pbc get returns one list per frame; -now is the current frame)
set cell_vec [pbc get -now -molid $topid]
set cell_size [lrange [lindex $cell_vec 0] 0 2]
set cell_shape [lrange [lindex $cell_vec 0] 3 5]
if { [lindex $cell_size 0] == 0 } {
    puts "WARNING: no box information found in the coordinates; set the &CELL of forceeval_subsys.cp2k.inc by hand!"
}
lappend subsys_list "&CELL\n  ABC \[angstrom\] [join $cell_size "  "]     !size of cell\n  ALPHA_BETA_GAMMA [join $cell_shape " "]      !shape of cell\n&END CELL\n"

#--> set QM layer
set qm_layer_sel [atomselect $topid $qm_layer frame first]
puts "CONTROL STRING: [$qm_layer_sel num] atoms were selected!"
if { [$qm_layer_sel num] == 0 } {
    puts "ERROR: the QM selection is empty."
    return
}

mol selection "noh $qm_layer"
mol representation DynamicBonds 1.7 0.1
mol material AOChalky
mol addrep $topid

mol selection "$qm_layer"
mol representation DynamicBonds 1.3 0.1
mol material AOChalky
mol addrep $topid

set qm_resi [lsort -unique -integer [$qm_layer_sel get resid]]
mol selection "resid $qm_resi and not ($qm_layer)"
mol representation Lines 3
mol material AOChalky
mol addrep $topid

#--> set INCLUDE file for CP2K QMMM section
set qmmm_list [list]

#--> set QM cell (extent of the QM layer plus the padding on each side)
lappend qmmm_list "! CP2K > FORCE_EVAL > QMMM > CELL"
set qmcell_vec [measure minmax $qm_layer_sel]
set pad2 [expr 2.0*$qm_pad]
set qmcell_size [vecadd [vecsub [lindex $qmcell_vec 1] [lindex $qmcell_vec 0]] [list $pad2 $pad2 $pad2]]
set qmcell_shape [list 90.0 90.0 90.0]
lappend qmmm_list "&CELL\n  ABC \[angstrom\] [join $qmcell_size "  "]     !size of cell\n  ALPHA_BETA_GAMMA [join $qmcell_shape " "]      !shape of cell\n&END CELL\n"

#--> set QM atoms (one &QM_KIND per element, guessed from the mass)
lappend qmmm_list "! CP2K > FORCE_EVAL > QMMM > QM_KIND"
set qm_masses [lsort -unique -real [$qm_layer_sel get mass]]
foreach mass $qm_masses {
    set kind_sel [atomselect $topid "mass $mass and ($qm_layer)" frame first]
    set qm_kind [lsort -unique [$kind_sel get element]]
    if { [llength $qm_kind] > 1 } {
        set qm_kind [lsearch -inline -all -not -exact $qm_kind X]
    }
    mol selection "mass $mass and ($qm_layer)"
    mol representation VDW 0.2
    mol material AOChalky
    mol addrep $topid
    lappend qmmm_list "&QM_KIND $qm_kind   !QM_KIND guess from element mass $mass\n  MM_INDEX [$kind_sel get serial]\n&END QM_KIND"
}

#--> set link atoms (every bond between a heavy QM atom and an MM atom)
lappend qmmm_list "! CP2K > FORCE_EVAL > QMMM > LINK"
set qm_index [$qm_layer_sel get index]
set heavy_qm_index [[atomselect $topid "noh $qm_layer" frame first] get index]
set nlink 0
foreach id $heavy_qm_index {
    foreach conect [lindex [[atomselect $topid "index $id"] getbonds] 0] {
        if { [lsearch -exact $qm_index $conect] == -1 } {
            incr nlink
            mol selection "index $id $conect"
            mol representation CPK 0.3 0.2
            mol material AOChalky
            mol addrep $topid
            label add Bonds $topid/$id $topid/$conect
            set lnk_bond [measure bond [list $id $conect] molid $topid]
            #rough guess for ALPHA_IMOMM calculated from the ratio between the cut LNK_BND and the CT-HC bond, applying a penalty for polarity on the CT-HC bond based on the difference between the CT-CT bond and the LNK_BND - check LNK_BND during minimization to validate
            set alpha_imomm [expr $lnk_bond/(1.090-(1.526-$lnk_bond))]
            lappend qmmm_list "&LINK\n  ALPHA_IMOMM [format "%.2f" $alpha_imomm]     ![join [[atomselect $topid "index $id $conect"] get type] "-"] bond (default is ALPHA_IMOMM 1.38 for C-C bond)\n  LINK_TYPE IMOMM\n  MM_INDEX  [expr $conect+1]\n  QM_INDEX  [expr $id+1]\n  QMMM_SCALE_FACTOR [format "%.1f" $qmmm_sc]\n&END LINK"
        }
    }
}
puts "CONTROL STRING: $nlink QM/MM bonds were cut (one &LINK section each)!"

#--> write INCLUDE file for CP2K SUBSYS section
set cp2k_inc [open forceeval_subsys.cp2k.inc $out_mode]
puts $cp2k_inc [join $subsys_list "\n"]
close $cp2k_inc

#--> write INCLUDE file for CP2K QMMM section
set cp2k_inc [open forceeval_qmmm.cp2k.inc $out_mode]
puts $cp2k_inc [join $qmmm_list "\n"]
close $cp2k_inc

#--> correct atom charge with parmed
set parmed_list [list]
set molname [molinfo $topid get name]
lappend parmed_list "parm $molname\n"
puts "\n3) DOUBLE-CHECK THAT INTEGER CHARGE OF SELECTION MATCHES ITS PARCELS!"
set qm_q [expr [join [$qm_layer_sel get charge] +]]
puts "CONTROL STRING: Charge of selection is $qm_q!"

foreach resi $qm_resi {
    set qmresi_q [expr [join [[atomselect $topid "resid $resi and ($qm_layer)" frame first] get charge] +]]
    set mm_atom [[atomselect $topid "resid $resi and not ($qm_layer)" frame first] get serial]
    if { [llength $mm_atom] > 0 } {
        #difference between atomic charges in QM residue and expected integer charge balance is distributed among residue MM atoms - make sure charge balance makes sense
        set dq [expr ($qmresi_q - round($qmresi_q))/[llength $mm_atom]]
        puts "CONTROL STRING: Charge of the selection in resid $resi is $qmresi_q (applying correction $dq to make it [expr round($qmresi_q)])!"

        foreach atom $mm_atom {
            set q [[atomselect $topid "serial $atom" frame first] get charge]
            lappend parmed_list "change charge @$atom [expr $q + $dq]"
        }
    }
}
#correct Lennard-Jones parameters for hydrogen in OH groups (OH,oh) and waters (HW,hw) - double-check for additional parameterized species
foreach t $lj_types {
    lappend parmed_list "changeLJSingleType @%$t 0.3019 0.047"
}
lappend parmed_list "\noutparm [file rootname $molname]_ee.prmtop"

#--> write SCRIPT for PARMED charge correction
set parmed_file [open qmmm-ee.parmed.in $out_mode]
puts $parmed_file [join $parmed_list "\n"]
close $parmed_file

puts "\nFiles written: forceeval_subsys.cp2k.inc forceeval_qmmm.cp2k.inc qmmm-ee.parmed.in"
puts "Next: parmed -i qmmm-ee.parmed.in   (writes [file rootname $molname]_ee.prmtop)"
if { $quit_at_end } { quit }
