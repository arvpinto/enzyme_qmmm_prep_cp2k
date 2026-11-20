#!/usr/bin/env bash
# mol2_qmsel.sh
#
# Converts a MOL2 file of the QM layer (e.g. saved from PyMOL with "save hl.mol2, sele")
# into the equivalent VMD, PyMOL and Amber mask selections, residue by residue.
# The MOL2 residue labels must be RESNAME+RESID (e.g. ASP196, 4GA500), as PyMOL writes them,
# and the residue numbers must be those of the Amber topology.
#
# Usage:
#   mol2_qmsel.sh hl_65atoms.mol2
#   mol2_qmsel.sh -o qm_sel hl_65atoms.mol2     # also writes qm_sel.vmd, qm_sel.pml and qm_sel.mask
#
# Options:
#   -o PREFIX   write each selection to PREFIX.vmd, PREFIX.pml and PREFIX.mask (one line each)
#   -h          print this help

usage() { sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; }

PREFIX=""
while getopts ":o:h" opt; do
    case $opt in
        o) PREFIX=$OPTARG ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done
shift $((OPTIND-1))

HL=$1
if [[ -z "$HL" || ! -f "$HL" ]]; then
    echo "ERROR: MOL2 file not given or not found." >&2
    usage
    exit 1
fi

#RESLIST=`sed -n '/SUBSTRUCTURE/,$p' $HL | grep -oE '^[0-9]+\s*\S+\s*' | cut -d$'\t' -f2`
RESLIST=`sed -n '/SUBSTRUCTURE/,$p' $HL | cut -d$'\t' -f2 | tail -n+2`
HLIST=`sed -n '/ATOM/,/BOND/p' $HL`

for RES in ${RESLIST[@]}; do
    NAME=( $(echo "$HLIST" | grep -w $RES | sed -E 's/^[0-9]+\s*//' | grep -oE '^\S+') )
    NAME=$(echo ${NAME[@]})
    RESN=$(echo $RES | cut -c1-3)
    RESI=$(echo $RES | cut -c4-)
    VMD_SEL=$VMD_SEL"`echo "(name ${NAME} and resname ${RESN} and resid ${RESI})"`"
    PML_SEL=$PML_SEL"`echo "(name ${NAME// /+} & resn ${RESN} & resi ${RESI})"`"
    AMB_SEL=$AMB_SEL"`echo "(@${NAME// /,}&:${RESN}&:${RESI})"`"
done
NATOM=$(echo "$HLIST" | grep -cE '^\s*[0-9]+\s')
echo -e "RESIDUE LIST ($NATOM atoms):\n$RESLIST"
echo -e "\nVMD selection:\n${VMD_SEL//)(/) or (}"
echo -e "\nPYMOL selection:\n${PML_SEL//)(/) | (}"
echo -e "\nAMBER MASK selection:\n${AMB_SEL//)(/)|(}"

if [[ -n "$PREFIX" ]]; then
    echo "${VMD_SEL//)(/) or (}" > $PREFIX.vmd
    echo "${PML_SEL//)(/) | (}" > $PREFIX.pml
    echo "${AMB_SEL//)(/)|(}" > $PREFIX.mask
    echo -e "\nSelections written to $PREFIX.vmd, $PREFIX.pml and $PREFIX.mask"
fi
