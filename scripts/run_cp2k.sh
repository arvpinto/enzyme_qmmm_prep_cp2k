#!/usr/bin/env bash
# run_cp2k.sh
#
# Runs the CP2K stages of the tutorial in order. Each stage runs in its own directory (as in the
# repository) and gets the files it reads copied there first: the topologies and include files from
# system_prep/, and the restart and wavefunction files from the directory of the previous stage.
# After each stage the script checks that CP2K ended normally.
#
#   stage        directory      input              reads (besides the input)
#   mm_em        mm_em/         mm_em.inp          amylase.prmtop, frame_xray.rst7, forceeval_subsys.cp2k.inc
#   qmmm_em_me   qmmm_em/       qmmm_em_me.inp     amylase.prmtop, forceeval_qmmm.cp2k.inc, MM_EM-1.restart
#   qmmm_em_ee   qmmm_em/       qmmm_em_ee.inp     amylase_ee.prmtop, forceeval_qmmm.cp2k.inc, QMMM_EM_ME-1.restart/-RESTART.wfn
#   qmmm_md_eq   qmmm_md_eq/    qmmm_md_eq.inp     amylase_ee.prmtop, forceeval_qmmm.cp2k.inc, QMMM_EM_EE-1.restart/-RESTART.wfn
#   qmmm_md_smd  qmmm_md_smd/   qmmm_md_smd.inp    amylase_ee.prmtop, forceeval_qmmm.cp2k.inc, QMMM_MD_EQ-1.restart(.bin)/-RESTART.wfn
#
# The repository already holds the outputs of every stage. To run them again without touching it,
# give a work directory with -w: the stage directories are created there, and files produced in the
# work directory are preferred to those of the repository.
#
# Usage:
#   ./scripts/run_cp2k.sh -h
#   ./scripts/run_cp2k.sh -w ~/scratch/amylase -s mm_em -n 64          # MM minimization in a work directory
#   ./scripts/run_cp2k.sh -w ~/scratch/amylase -n 128 -l srun          # all stages, inside a Slurm job
#   ./scripts/run_cp2k.sh -w ~/scratch/amylase -s qmmm_md_smd -c       # only copy the files of a stage
#   ./scripts/run_cp2k.sh -d                                           # print what would be done
#
# Options:
#   -s "STAGES"   stages to run, in this order (default: "mm_em qmmm_em_me qmmm_em_ee qmmm_md_eq qmmm_md_smd")
#   -w DIR        work directory (default: the repository itself)
#   -r DIR        repository root (default: the parent directory of this script)
#   -n N          number of MPI tasks (default: 4)
#   -t N          OpenMP threads per MPI task, exported as OMP_NUM_THREADS (default: 1)
#   -e EXE        CP2K executable (default: cp2k.psmp)
#   -l LAUNCHER   MPI launcher: mpirun (default), srun, or none to run EXE directly
#   -c            only copy the files each stage needs, do not run CP2K
#   -f            run a stage even if its output already finished (the old output is renamed .bak)
#   -d            dry run: report the files and print the commands without copying or running
#   -h            print this help

usage() { sed -n '2,37p' "$0" | sed 's/^# \{0,1\}//'; }

STAGES="mm_em qmmm_em_me qmmm_em_ee qmmm_md_eq qmmm_md_smd"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK=""
NP=4
NT=1
EXE=cp2k.psmp
LAUNCHER=mpirun
COPY_ONLY=0
FORCE=0
DRY=0
while getopts ":s:w:r:n:t:e:l:cfdh" opt; do
    case $opt in
        s) STAGES=$OPTARG ;;
        w) WORK=$OPTARG ;;
        r) ROOT=$OPTARG ;;
        n) NP=$OPTARG ;;
        t) NT=$OPTARG ;;
        e) EXE=$OPTARG ;;
        l) LAUNCHER=$OPTARG ;;
        c) COPY_ONLY=1 ;;
        f) FORCE=1 ;;
        d) DRY=1 ;;
        h) usage; exit 0 ;;
        :) echo "ERROR: option -$OPTARG needs a value" >&2; usage; exit 1 ;;
        *) echo "ERROR: unknown option -$OPTARG" >&2; usage; exit 1 ;;
    esac
done
[[ -z "$WORK" ]] && WORK=$ROOT
if [[ ! -d "$ROOT/system_prep" ]]; then
    echo "ERROR: $ROOT does not look like the repository (no system_prep/ directory); use -r" >&2
    exit 1
fi
mkdir -p "$WORK" 2>/dev/null || true
WORK="$(cd "$WORK" && pwd)"
export OMP_NUM_THREADS=$NT

# stage -> directory, input, and the files it reads ("dir/file", relative to the repository root)
declare -A DIR=( [mm_em]=mm_em [qmmm_em_me]=qmmm_em [qmmm_em_ee]=qmmm_em [qmmm_md_eq]=qmmm_md_eq [qmmm_md_smd]=qmmm_md_smd )
declare -A NEEDS=(
    [mm_em]="system_prep/amylase.prmtop system_prep/frame_xray.rst7 system_prep/forceeval_subsys.cp2k.inc"
    [qmmm_em_me]="system_prep/amylase.prmtop system_prep/forceeval_qmmm.cp2k.inc mm_em/MM_EM-1.restart"
    [qmmm_em_ee]="system_prep/amylase_ee.prmtop system_prep/forceeval_qmmm.cp2k.inc qmmm_em/QMMM_EM_ME-1.restart qmmm_em/QMMM_EM_ME-RESTART.wfn"
    [qmmm_md_eq]="system_prep/amylase_ee.prmtop system_prep/forceeval_qmmm.cp2k.inc qmmm_em/QMMM_EM_EE-1.restart qmmm_em/QMMM_EM_EE-RESTART.wfn"
    [qmmm_md_smd]="system_prep/amylase_ee.prmtop system_prep/forceeval_qmmm.cp2k.inc qmmm_md_eq/QMMM_MD_EQ-1.restart qmmm_md_eq/QMMM_MD_EQ-1.restart.bin qmmm_md_eq/QMMM_MD_EQ-RESTART.wfn" )

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# a file produced in the work directory wins over the copy in the repository
source_of() {
    if [[ -f "$WORK/$1" ]]; then echo "$WORK/$1"
    elif [[ -f "$ROOT/$1" ]]; then echo "$ROOT/$1"
    fi
}

if [[ $DRY -eq 0 && $COPY_ONLY -eq 0 ]] && ! command -v "$EXE" > /dev/null; then
    echo "ERROR: CP2K executable $EXE not found in PATH (use -e)" >&2
    exit 1
fi
case "$LAUNCHER" in
    none) RUN="$EXE" ;;
    srun) RUN="srun -n $NP $EXE" ;;
    *)    RUN="$LAUNCHER -np $NP $EXE" ;;
esac
log "repository: $ROOT"
log "work directory: $WORK"

for st in $STAGES; do
    d=${DIR[$st]}
    if [[ -z "$d" ]]; then
        echo "ERROR: unknown stage '$st' (valid: mm_em qmmm_em_me qmmm_em_ee qmmm_md_eq qmmm_md_smd)" >&2
        exit 1
    fi
    inp=$st.inp
    out=$st.out
    wd="$WORK/$d"
    [[ $DRY -eq 0 ]] && mkdir -p "$wd"

    # a stage that already finished is skipped before anything is copied
    if [[ -f "$wd/$out" && $FORCE -eq 0 && $DRY -eq 0 && $COPY_ONLY -eq 0 ]]; then
        if grep -q "PROGRAM ENDED AT" "$wd/$out"; then
            log "$st: $wd/$out already finished, skipping (use -f to run it again)"
            continue
        fi
        echo "ERROR: $wd/$out exists but did not finish; remove it or use -f" >&2
        exit 1
    fi

    # input file and the files the stage reads
    missing=""
    for f in "$d/$inp" ${NEEDS[$st]}; do
        src=$(source_of "$f")
        dst="$wd/$(basename "$f")"
        if [[ -z "$src" ]]; then
            missing="$missing $f"
        elif [[ "$src" != "$dst" ]]; then
            if [[ $DRY -eq 1 ]]; then
                log "$st: would copy $src -> $dst"
            else
                cp -p "$src" "$dst"
            fi
        fi
    done
    if [[ -n "$missing" ]]; then
        if [[ $DRY -eq 1 ]]; then
            log "$st: not available yet:$missing"
        else
            echo "ERROR: stage $st cannot start, missing:$missing" >&2
            exit 1
        fi
    fi
    if [[ $COPY_ONLY -eq 1 ]]; then
        log "$st: files copied to $wd"
        continue
    fi

    log "$st: cd $wd && $RUN -i $inp -o $out   (OMP_NUM_THREADS=$NT)"
    [[ $DRY -eq 1 ]] && continue
    [[ -f "$wd/$out" ]] && mv "$wd/$out" "$wd/$out.bak"
    ( cd "$wd" && $RUN -i "$inp" -o "$out" )
    if ! grep -q "PROGRAM ENDED AT" "$wd/$out"; then
        echo "ERROR: $st did not finish normally, check $wd/$out" >&2
        exit 1
    fi
    if grep -q "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" "$wd/$out"; then
        log "$st: WARNING, the optimization reached MAX_ITER without converging"
    fi
    log "$st: finished"
done
