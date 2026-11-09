#!/usr/bin/env vmd
# vmd_cp2k-repr.tcl
#
# Loads an Amber topology and CP2K trajectories (DCD) and draws the QM layer with dynamic bonds,
# so that bonds that form or break along a QM/MM MD run appear and disappear in the animation.
# Optionally labels distances (e.g. the atoms of a collective variable) to follow them in time.
#
# Before loading, the DCD headers are read to estimate how much memory the trajectory needs
# (VMD keeps ~12 bytes per atom per frame in RAM). If it does not fit in the available memory,
# frames are skipped automatically (stride) unless -step is given explicitly.
#
# Usage:
#   vmd -e vmd_cp2k-repr.tcl -args amylase_ee.prmtop QMMM_MD_SMD-pos-1.dcd
#       -> the script asks for the QM selection (VMD syntax) in the terminal
#   vmd -e vmd_cp2k-repr.tcl -args amylase_ee.prmtop QMMM_MD_SMD-pos-1.dcd \
#       -qminc forceeval_qmmm.cp2k.inc -bonds "7699:7700 3607:7699 3026:7700" -cartoon
#   vmd -e vmd_cp2k-repr.tcl -args amylase_ee.prmtop QMMM_MD_SMD-pos-1.dcd \
#       -qminc forceeval_qmmm.cp2k.inc -first 2000 -last 6000 -step 5
#
# Options (after -args; the .prmtop and the trajectory files can be given in any order):
#   <file>.prmtop        Amber topology (required)
#   <file>.dcd ...       trajectories, loaded in the order given (.xyz and .pdb files are also accepted)
#   -sel "<selection>"   QM layer in VMD selection syntax
#   -selfile <file>      file with the QM layer selection (VMD syntax)
#   -qminc <file>        take the QM layer from the &QM_KIND sections of a CP2K include/input file
#   -bonds "i:j ..."     pairs of atom serial numbers (1-based, as in CP2K) to label as distances
#                        (space- or comma-separated; quoting the list is recommended)
#   -first <n>           first frame to load from each trajectory (0-based, default 0)
#   -last <n>            last frame to load from each trajectory (default -1 = last frame)
#   -step <n>            load every n-th frame (default: 1, or larger if the trajectory
#                        would not fit in memory)
#   -maxmem <GB>         memory allowed for the coordinates (default: 60% of the available RAM)
#   -nocheck             skip the memory check and load with the given -step (default 1)
#   -cartoon             also show the protein as a cartoon (created hidden otherwise)
#   -h                   print this help

proc repr_usage {} {
    puts "Usage: vmd -e vmd_cp2k-repr.tcl -args <prmtop> <traj.dcd> \[...\] \[-sel \"<sel>\" | -selfile <file> | -qminc <file>\]"
    puts "       \[-bonds \"i:j ...\"\] \[-first n\] \[-last n\] \[-step n\] \[-maxmem GB\] \[-nocheck\] \[-cartoon\]"
}

# Available memory in bytes (Linux /proc/meminfo), or -1 if unknown (e.g. macOS)
proc mem_available_bytes {} {
    if { ![file readable /proc/meminfo] } { return -1 }
    set fh [open /proc/meminfo r]; set txt [read $fh]; close $fh
    if { [regexp -line {^MemAvailable:\s+(\d+)\s+kB} $txt -> kb] } { return [expr {wide($kb) * 1024}] }
    if { [regexp -line {^MemFree:\s+(\d+)\s+kB} $txt -> kb] }      { return [expr {wide($kb) * 1024}] }
    return -1
}

# Read a CHARMM/NAMD/CP2K DCD header; returns {natoms nframes} or {} if not a readable DCD
proc dcd_info { fname } {
    if { [catch {open $fname r} fh] } { return {} }
    fconfigure $fh -translation binary
    set hdr [read $fh 92]
    if { [string length $hdr] < 92 } { close $fh; return {} }
    set e ""
    if { [binary scan $hdr i m] && $m == 84 } { set e i } elseif { [binary scan $hdr I m] && $m == 84 } { set e I }
    if { $e eq "" || [string range $hdr 4 7] ne "CORD" } { close $fh; return {} }
    binary scan $hdr x8${e}20 icntrl
    set has_cell [expr {[lindex $icntrl 10] != 0}]
    set has_4d   [expr {[lindex $icntrl 11] != 0}]
    # title record
    binary scan [read $fh 4] $e tlen
    seek $fh [expr {$tlen + 4}] current
    # natoms record
    binary scan [read $fh 12] ${e}3 natrec
    set natoms [lindex $natrec 1]
    set hdrsize [tell $fh]
    close $fh
    if { $natoms <= 0 } { return {} }
    set framesize [expr {wide(12) * $natoms + 24}]
    if { $has_cell } { set framesize [expr {$framesize + 56}] }
    if { $has_4d }   { set framesize [expr {$framesize + 4 * $natoms + 8}] }
    set nframes [expr {([file size $fname] - $hdrsize) / $framesize}]
    return [list $natoms $nframes]
}

# Number of frames VMD will load from a file with nframes frames for first/last/step
proc frames_in_range { nframes first last step } {
    if { $last < 0 || $last >= $nframes } { set last [expr {$nframes - 1}] }
    if { $first > $last } { return 0 }
    return [expr {($last - $first) / $step + 1}]
}

set topol ""
set trajs [list]
set qm_str ""
set bond_pairs [list]
set first 0
set last -1
set step ""
set maxmem_gb ""
set memcheck 1
set show_cartoon 0
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
                set l [string trim [lindex [split $line "!#"] 0]]
                if { [string match -nocase "&QM_KIND*" $l] } { set in_kind 1; continue }
                if { [string match -nocase "&END*" $l] } { set in_kind 0; continue }
                if { $in_kind } {
                    regsub -nocase {^MM_INDEX} $l "" l
                    regsub -all {\\} $l " " l
                    foreach s $l { if { [string is integer -strict $s] } { lappend serials $s } }
                }
            }
            if { [llength $serials] == 0 } {
                puts "ERROR: no MM_INDEX entries found in the &QM_KIND sections of [lindex $argv $i]."
                return
            }
            set qm_str "serial [lsort -unique -integer $serials]"
        }
        -bonds   {
            # accept "i:j k:l" as one argument, and also i:j k:l as separate arguments
            # (in case the quoted string was split on the way into VMD's argv)
            incr i
            foreach p [split [string map {, " "} [lindex $argv $i]]] { if { $p ne "" } { lappend bond_pairs $p } }
            while { $i + 1 < [llength $argv] && [regexp {^\d+[:-]\d+$} [lindex $argv [expr {$i + 1}]]] } {
                incr i
                lappend bond_pairs [lindex $argv $i]
            }
        }
        -first   { incr i; set first [lindex $argv $i] }
        -last    { incr i; set last [lindex $argv $i] }
        -step    { incr i; set step [lindex $argv $i] }
        -maxmem  { incr i; set maxmem_gb [lindex $argv $i] }
        -nocheck { set memcheck 0 }
        -cartoon { set show_cartoon 1 }
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
foreach {name val} [list -first $first -last $last] {
    if { ![string is integer -strict $val] } { puts "ERROR: $name must be an integer."; return }
}
if { $step ne "" && (![string is integer -strict $step] || $step < 1) } {
    puts "ERROR: -step must be a positive integer."; return
}
if { $maxmem_gb ne "" && (![string is double -strict $maxmem_gb] || $maxmem_gb <= 0) } {
    puts "ERROR: -maxmem must be a positive number (GB)."; return
}

display rendermode GLSL
puts "\n1) MAKE SURE YOUR SYSTEM IS THE TOP MOLECULE!"
mol new $topol type parm7 first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
set topid [molinfo top get id]
set natoms [molinfo $topid get numatoms]
puts "CONTROL STRING: TOP MOLECULE IS $topid ($natoms atoms)!"

if { $qm_str eq "" } {
    puts "\n2) INTRODUCE SELECTION FOR QM LAYER IN VMD SELECTION FORMAT:"
    set qm_str [gets stdin]
} else {
    puts "\n2) QM LAYER SELECTION:\n$qm_str"
}
set qm_sel [atomselect $topid $qm_str frame first]
if { [$qm_sel num] == 0 } {
    puts "ERROR: the QM selection is empty, check it against the topology."
    $qm_sel delete
    return
}
set qm_resi [lsort -unique -integer [$qm_sel get resid]]
set qm_q [expr [join [$qm_sel get charge] +]]
set qm_index [$qm_sel get index]
set heavy_sel [atomselect $topid "noh ($qm_str)" frame first]
set heavyqm_index [$heavy_sel get index]
set heavyqm_bonds [$heavy_sel getbonds]
$heavy_sel delete
puts "CONTROL STRING: [$qm_sel num] atoms were selected, with total charge of selection is [format %.4f $qm_q]!"
$qm_sel delete

puts "\n3) LOAD TRAJECTORIES ON SYSTEM TOPOLOGY:"
# --- memory check: estimate the RAM needed for the coordinates before loading
set bytes_per_frame [expr {wide(12) * $natoms * 1.1}]
set nframes_range 0
set unknown_size 0
foreach trj $trajs {
    if { ![file exists $trj] } { puts "ERROR: trajectory $trj not found."; return }
    set info [dcd_info $trj]
    if { [llength $info] == 0 } { set unknown_size 1; continue }
    lassign $info dcd_natoms dcd_nframes
    if { $dcd_natoms != $natoms } {
        puts "ERROR: $trj has $dcd_natoms atoms but $topol has $natoms atoms."
        return
    }
    set n [frames_in_range $dcd_nframes $first $last 1]
    puts "CONTROL STRING: $trj contains $dcd_nframes frames ($n in the requested range)."
    incr nframes_range $n
}
set avail [mem_available_bytes]
if { $maxmem_gb ne "" } {
    set limit [expr {$maxmem_gb * 1024.0 * 1024 * 1024}]
} elseif { $avail > 0 } {
    set limit [expr {0.6 * $avail}]
} else {
    set limit -1
}
set need_gb [expr {$nframes_range * $bytes_per_frame / 1024.0**3}]
if { $memcheck && $limit > 0 && $nframes_range > 0 } {
    set limit_gb [expr {$limit / 1024.0**3}]
    set min_step [expr {int(ceil($nframes_range * $bytes_per_frame / $limit))}]
    if { $min_step < 1 } { set min_step 1 }
    puts [format "CONTROL STRING: all frames would need %.2f GB; memory limit is %.2f GB." $need_gb $limit_gb]
    if { $step eq "" } {
        set step $min_step
        if { $step > 1 } {
            puts "WARNING: not enough memory for every frame, loading every ${step}th frame (use -first/-last to zoom into a window at full resolution, or -maxmem to change the limit)."
        }
    } elseif { $step < $min_step } {
        puts "WARNING: with -step $step the trajectory needs about [format %.2f [expr {$need_gb / $step}]] GB; VMD may run out of memory. Use -step $min_step or larger."
    }
} else {
    if { $step eq "" } { set step 1 }
    if { $unknown_size || $limit <= 0 } {
        puts "CONTROL STRING: memory check skipped (non-DCD trajectory or unknown available memory)."
    }
}
foreach trj $trajs {
    set ext [string tolower [string trimleft [file extension $trj] "."]]
    if { $ext eq "" } { set ext dcd }
    mol addfile $trj type $ext first $first last $last step $step filebonds 1 autobonds 1 waitfor all molid $topid
}
puts "CONTROL STRING: Trajectories $trajs were loaded on $topol ([molinfo $topid get numframes] frames, first $first, last $last, step $step)!"

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

#protein cartoon, hidden unless -cartoon is given (switch it on in Graphics > Representations)
mol selection "protein and not resid $qm_resi"
mol representation NewCartoon 0.3 10.0 4.1 0
mol material AOChalky
mol color ColorID 6
mol addrep $topid
if { !$show_cartoon } { mol showrep $topid [expr [molinfo $topid get numreps] - 1] 0 }

set link_index [list]
foreach id $heavyqm_index conects $heavyqm_bonds {
    foreach conect $conects {
        if { [lsearch -exact $qm_index $conect] == -1 } {
            lappend link_index "$id" "$conect"
            if { [catch {label add Bonds $topid/$id $topid/$conect} err] } {
                puts "WARNING: could not label the QM/MM boundary bond $id-$conect: $err"
            }
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
    if { ![regexp {^(\d+)[:-](\d+)$} [string trim $pair] -> s1 s2] } {
        puts "WARNING: -bonds entry \"$pair\" is not of the form i:j, skipped."
        continue
    }
    if { $s1 < 1 || $s2 < 1 || $s1 > $natoms || $s2 > $natoms } {
        puts "WARNING: -bonds pair $s1:$s2 is outside the serial range 1-$natoms, skipped."
        continue
    }
    if { $s1 == $s2 } {
        puts "WARNING: -bonds pair $s1:$s2 uses the same atom twice, skipped."
        continue
    }
    set a1 [expr {$s1 - 1}]
    set a2 [expr {$s2 - 1}]
    if { [catch {label add Bonds $topid/$a1 $topid/$a2} err] } {
        puts "WARNING: could not add the label for serials $s1:$s2 (indices $a1:$a2): $err"
        continue
    }
    set names [list]
    foreach a [list $a1 $a2] {
        set s [atomselect $topid "index $a"]
        lappend names "[$s get resname][$s get resid]:[$s get name]"
        $s delete
    }
    puts "CONTROL STRING: label for the distance between serials $s1 and $s2 ([join $names { - }]) (Graphics > Labels > Bonds > Graph)"
}

animate goto 0
display resetview
puts "CONTROL STRING: End of molecule representations section!"
