#!/usr/bin/env vmd
# vmd_cp2k-repr.tcl
#
# Loads an Amber topology and CP2K trajectories (DCD) and draws the QM layer with dynamic bonds,
# so that bonds that form or break along a QM/MM MD run appear and disappear in the animation.
# Optionally labels distances (e.g. the atoms of a collective variable) to follow them in time.
#
# Usage:
#   vmd -e vmd_cp2k-repr.tcl -args amylase_ee.prmtop SMD-NVT-EE-GPW-pos-1.dcd
#       -> the script asks for the QM selection (VMD syntax) in the terminal
#   vmd -e vmd_cp2k-repr.tcl -args amylase_ee.prmtop REACT-NVT-EE-GPW-pos-1.dcd SMD-NVT-EE-GPW-pos-1.dcd \
#       -qminc forceeval_qmmm.cp2k.inc -bonds "7699:7700 3607:7699 3026:7700"
#
# Options (after -args; the .prmtop and the trajectory files can be given in any order):
#   <file>.prmtop        Amber topology (required)
#   <file>.dcd ...       trajectories, loaded in the order given (.xyz and .pdb files are also accepted)
#   -sel "<selection>"   QM layer in VMD selection syntax
#   -selfile <file>      file with the QM layer selection (VMD syntax)
#   -qminc <file>        take the QM layer from the &QM_KIND sections of a CP2K include/input file
#   -bonds "i:j ..."     pairs of atom serial numbers (1-based, as in CP2K) to label as distances
#   -step <n>            load every n-th frame (default 1)
#   -h                   print this help

proc repr_usage {} {
    puts "Usage: vmd -e vmd_cp2k-repr.tcl -args <prmtop> <traj.dcd> \[...\] \[-sel \"<sel>\" | -selfile <file> | -qminc <file>\] \[-bonds \"i:j ...\"\] \[-step n\]"
}

set topol ""
set trajs [list]
set qm_str ""
set bond_pairs [list]
set step 1
for {set i 0} {$i < [llength $argv]} {incr i} {
    set opt [lindex $argv $i]
    switch -glob -- $opt {
        -sel     { incr i; set qm_str [lindex $argv $i] }
        -selfile { incr i; set fh [open [lindex $argv $i] r]; set qm_str [string trim [join [split [read $fh] "\n"] " "]]; close $fh }
        -qminc   {
            incr i
            set fh [open [lindex $argv $i] r]; set txt [read $fh]; close $fh
            set serials [list]
            set in_kind 0
            foreach line [split $txt "\n"] {
                set l [string trim [lindex [split $line "!"] 0]]
                if { [string match -nocase "&QM_KIND*" $l] } { set in_kind 1; continue }
                if { [string match -nocase "&END*" $l] } { set in_kind 0; continue }
                if { $in_kind } {
                    regsub -nocase {^MM_INDEX} $l "" l
                    regsub -all {\\} $l " " l
                    foreach s $l { if { [string is integer -strict $s] } { lappend serials $s } }
                }
            }
            set qm_str "serial $serials"
        }
        -bonds   { incr i; set bond_pairs [lindex $argv $i] }
        -step    { incr i; set step [lindex $argv $i] }
        -h       { repr_usage; return }
        *.prmtop - *.parm7 { set topol $opt }
        -*       { puts "Unknown option $opt"; repr_usage; return }
        default  { lappend trajs $opt }
    }
}
if { $topol eq "" } {
    puts "ERROR: no Amber topology (.prmtop) given."
    repr_usage
    return
}

display rendermode GLSL
puts "\n1) MAKE SURE YOUR SYSTEM IS THE TOP MOLECULE!"
mol new $topol type parm7 first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
set topid [molinfo top get id]
puts "CONTROL STRING: TOP MOLECULE IS $topid!"

if { $qm_str eq "" } {
    puts "\n2) INTRODUCE SELECTION FOR QM LAYER IN VMD SELECTION FORMAT:"
    set qm_str [gets stdin]
} else {
    puts "\n2) QM LAYER SELECTION:\n$qm_str"
}
set qm_sel [atomselect $topid $qm_str frame first]
set qm_resi [lsort -unique -integer [$qm_sel get resid]]
set qm_q [expr [join [$qm_sel get charge] +]]
set qm_index [$qm_sel get index]
set heavyqm_index [[atomselect $topid "noh ($qm_str)" frame first] get index]
puts "CONTROL STRING: [$qm_sel num] atoms were selected, with total charge of selection is $qm_q!"

puts "\n3) LOAD TRAJECTORIES ON SYSTEM TOPOLOGY:"
foreach trj $trajs {
    set ext [string tolower [string trimleft [file extension $trj] "."]]
    if { $ext eq "" } { set ext dcd }
    mol addfile $trj type $ext first 0 last -1 step $step filebonds 1 autobonds 1 waitfor all molid $topid
}
puts "CONTROL STRING: Trajectories $trajs were loaded on $topol ([molinfo $topid get numframes] frames)!"

puts "\n4) CREATE VMD MOLECULAR REPRESENTATIONS:"
mol delrep 0 $topid

mol selection "noh ($qm_str)"
mol representation DynamicBonds 1.8 0.1
mol material AOChalky
mol addrep $topid

mol selection "$qm_str"
mol representation DynamicBonds 1.3 0.1
mol material AOChalky
mol addrep $topid

mol selection "element H P S and ($qm_str)"
mol representation DynamicBonds 1.5 0.1
mol material AOChalky
mol addrep $topid

mol selection "element C P S and ($qm_str)"
mol representation DynamicBonds 2.1 0.1
mol material AOChalky
mol addrep $topid

mol selection "element P S and ($qm_str)"
mol representation DynamicBonds 2.4 0.1
mol material AOChalky
mol addrep $topid

mol selection "resid $qm_resi and not ($qm_str)"
mol representation Lines 3
mol material AOChalky
mol addrep $topid

mol selection "$qm_str"
mol representation VDW 0.2
mol material AOChalky
mol addrep $topid

#protein cartoon, created hidden (switch it on in Graphics > Representations)
mol selection "protein and not resid $qm_resi"
mol representation NewCartoon 0.3 10.0 4.1 0
mol material AOChalky
mol color ColorID 6
mol addrep $topid
mol showrep $topid [expr [molinfo $topid get numreps] - 1] 0

set link_index [list]
foreach id $heavyqm_index {
    foreach conect [lindex [[atomselect $topid "index $id"] getbonds] 0] {
        if { [lsearch -exact $qm_index $conect] == -1 } {
            lappend link_index "$id" "$conect"
            label add Bonds $topid/$id $topid/$conect
        }
    }
}
if { [llength $link_index] > 0 } {
    mol selection "index $link_index"
    mol representation CPK 0.3 0.2
    mol material AOChalky
    mol color Name
    mol addrep $topid
}

#--> distances to follow along the trajectory (serial numbers, 1-based)
foreach pair $bond_pairs {
    lassign [split $pair ":"] s1 s2
    label add Bonds $topid/[expr $s1-1] $topid/[expr $s2-1]
    puts "CONTROL STRING: label for the distance between serials $s1 and $s2 (Graphics > Labels > Bonds > Graph)"
}

animate goto 0
display resetview
puts "CONTROL STRING: End of molecule representations section!"
