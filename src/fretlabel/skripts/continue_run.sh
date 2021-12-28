#!/bin/bash

usage() { echo "Continue a stopped run 
Usage: continue_run.sh -s <prev. run input file (.tpr)> -c <prev. checkpoint file>

or extend a finished run
Usage: continue_run.sh -s <prev. run input file (.tpr)> -c <prev. checkpoint file> -e <timetoextend (ps)> -o <extension run input file (.tpr)>" 1>&2; exit 1; }
invalidOpt() { echo "Invalid option: -$OPTARG" 1>&2; exit 1; }
missingArg() { echo "Option -$OPTARG requires an argument" 1>&2; exit 1; }
cleanup() { if ls -f $1/\#* 1> /dev/null 2>&1 ; then rm $1/\#* ; fi ; }

while getopts ":s:c:e:o:h" opt; do
    case $opt in
        s) 
            inputFile=$OPTARG
            ;;
        c)
            checkptFile=$OPTARG
            ;;
        e)
            timetoextend=$OPTARG
            ;;
        o)
            extInputFile=$OPTARG
            ;;
        h)
            usage
            ;;
        \?)
            invalidOpt
            ;;
        :)
            missingArg
            ;;
        *)
            usage
            ;;
    esac
done

# no cmd line arguments given
if [ -z "$inputFile" ] || [ -z "$checkptFile" ]; then
    usage
fi

inputName=`echo $inputFile | rev | cut -f1 -d"/" | rev | cut -f1 -d"."`


if [ -z "$timetoextend" ]; then
    # continue interrupted run
    gmx mdrun -v -s md0/"$inputName".tpr -cpi md0/"$inputName".cpt -x md0/"$inputName".xtc -e md0/"$inputName".edr -g md0/"$inputName".log -append -cpo md0/"$inputName".cpt -c md0/"$inputName".gro || { echo "-> Error: gmx mdrun for md0 restart failed" ; cleanup md0 ; exit 1; }
else
    # extend finished run
    if [ -z "$extInputFile" ]; then
        usage
    fi
    gmx convert-tpr -s md0/"$inputName".tpr -extend "$timetoextend" -o md0/"$extInputFile" || { echo "-> Error: gmx convert-tpr for md0 extension failed" ; cleanup md0 ; exit 1; }
    gmx mdrun -v -s md0/"$extInputFile" -cpi md0/"$inputName".cpt -x md0/"$inputName".xtc -e md0/"$inputName".edr -g md0/"$inputName".log -cpo md0/"$inputName".cpt -c md0/"$inputName".gro || { echo "-> Error: gmx mdrun for md0 extension failed" ; cleanup md0 ; exit 1; }

fi 

cleanup md0
