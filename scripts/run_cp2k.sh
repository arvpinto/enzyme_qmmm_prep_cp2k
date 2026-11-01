#!/usr/bin/env bash
# run_cp2k.sh
#
# Runs the CP2K stages of the tutorial in order, checking before each stage that the files it
# reads (topology, include files, restart and wavefunction files of the previous stage) exist,
# and after each stage that CP2K ended normally.
#
#   em             em.inp                 MM minimization (FIST)
#   em-me          em-me_gpw.inp          QM/MM minimization, mechanical embedding
#   em-ee          em-ee_gpw.inp          QM/MM minimization, electrostatic embedding
#   react          react_nvt-ee_gpw.inp   QM/MM MD, NVT 300 K (reactant)
#   smd            smd-nvt-ee_gpw.inp     steered QM/MM MD along the CV
#
# Usage:
#   ./run_cp2k.sh                                   # all stages, mpirun -np 4 cp2k.psmp
#   ./run_cp2k.sh -n 32 -e cp2k.psmp -s "em-ee react smd"
#   ./run_cp2k.sh -l "srun" -n 128 -t 1             # inside a Slurm job, OMP_NUM_THREADS=1
#   ./run_cp2k.sh -d                                # print the commands only
#
# Options:
#   -s "STAGES"   stages to run, in this order (default: "em em-me em-ee react smd")
#   -n N          number of MPI tasks (default: 4)
#   -t N          OpenMP threads per MPI task, exported as OMP_NUM_THREADS (default: 1)
#   -e EXE        CP2K executable (default: cp2k.psmp)
#   -l LAUNCHER   MPI launcher (default: mpirun); "none" runs EXE without a launcher
#   -f            run a stage even if its output file already exists (it is renamed to .bak)
#   -d            dry run: check the files and print the commands without running them
#   -h            print this help

usage() { sed -n '2,31p' "$0" | sed 's/^# \{0,1\}//'; }

STAGES="em em-me em-ee react smd"
NP=4
NT=1
EXE=cp2k.psmp
LAUNCHER=mpirun
FORCE=0
DRY=0
while getopts ":s:n:t:e:l:fdh" opt; do
    case $opt in
        s) STAGES=$OPTARG ;;
        n) NP=$OPTARG ;;
        t) NT=$OPTARG ;;
        e) EXE=$OPTARG ;;
        l) LAUNCHER=$OPTARG ;;
        f) FORCE=1 ;;
        d) DRY=1 ;;
        h) usage; exit 0 ;;
        :) echo "ERROR: option -$OPTARG needs a value" >&2; usage; exit 1 ;;
        *) echo "ERROR: unknown option -$OPTARG" >&2; usage; exit 1 ;;
    esac
done
export OMP_NUM_THREADS=$NT

# stage -> input file, files it needs and file that marks a finished run
declare -A INP=( [em]=em.inp [em-me]=em-me_gpw.inp [em-ee]=em-ee_gpw.inp
                 [react]=react_nvt-ee_gpw.inp [smd]=smd-nvt-ee_gpw.inp )
declare -A NEEDS=(
    [em]="amylase.prmtop frame_xray.inpcrd forceeval_subsys.cp2k.inc"
    [em-me]="amylase.prmtop forceeval_qmmm.cp2k.inc EM-1.restart"
    [em-ee]="amylase_ee.prmtop forceeval_qmmm.cp2k.inc EM-ME-GPW-1.restart EM-ME-GPW-RESTART.wfn"
    [react]="amylase_ee.prmtop forceeval_qmmm.cp2k.inc EM-EE-GPW-1.restart EM-EE-GPW-RESTART.wfn"
    [smd]="amylase_ee.prmtop forceeval_qmmm.cp2k.inc REACT-NVT-EE-GPW-1.restart REACT-NVT-EE-GPW-1.restart.bin REACT-NVT-EE-GPW-RESTART.wfn" )

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

if [[ $DRY -eq 0 ]] && ! command -v $EXE > /dev/null; then
    echo "ERROR: CP2K executable $EXE not found in PATH (use -e)" >&2
    exit 1
fi
if [[ "$LAUNCHER" == "none" ]]; then
    RUN="$EXE"
elif [[ "$LAUNCHER" == "srun" ]]; then
    RUN="srun -n $NP $EXE"
else
    RUN="$LAUNCHER -np $NP $EXE"
fi

for st in $STAGES; do
    inp=${INP[$st]}
    if [[ -z "$inp" ]]; then
        echo "ERROR: unknown stage '$st' (valid: ${!INP[*]})" >&2
        exit 1
    fi
    out=${inp%.inp}.out
    missing=""
    for f in $inp ${NEEDS[$st]}; do
        [[ -f $f ]] || missing="$missing $f"
    done
    if [[ -n "$missing" ]]; then
        if [[ $DRY -eq 1 ]]; then
            log "$st: files not present yet:$missing"
        else
            echo "ERROR: stage $st cannot start, missing:$missing" >&2
            exit 1
        fi
    fi
    if [[ -f $out && $FORCE -eq 0 && $DRY -eq 0 ]]; then
        if grep -q "PROGRAM ENDED AT" $out; then
            log "$st: $out already finished, skipping (use -f to run it again)"
            continue
        fi
        echo "ERROR: $out exists but did not finish; remove it or use -f" >&2
        exit 1
    fi
    [[ -f $out && $DRY -eq 0 ]] && mv $out $out.bak
    log "$st: $RUN -i $inp -o $out   (OMP_NUM_THREADS=$NT)"
    [[ $DRY -eq 1 ]] && continue
    $RUN -i $inp -o $out
    if ! grep -q "PROGRAM ENDED AT" $out; then
        echo "ERROR: $st did not finish normally, check $out" >&2
        exit 1
    fi
    if [[ $st == em* ]] && grep -q "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" $out; then
        log "$st: WARNING, the optimization reached MAX_ITER without converging"
    fi
    log "$st: finished"
done
