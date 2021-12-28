#!/bin/bash

usage() { echo "Run a single replicate MD simulation
Usage: single_run.sh -c <structure file (.gro)> -d <mdp directory> -a <max. warnings for grompp> (optional) -p <plumed file (.dat)> -n <index file (.ndx, optional)>" 1>&2; exit 1; }
invalidOpt() { echo "Invalid option: -$OPTARG" 1>&2; exit 1; }
missingArg() { echo "Option -$OPTARG requires an argument" 1>&2; exit 1; }
cleanup() { if ls -f $1/\#* 1> /dev/null 2>&1 ; then rm $1/\#* ; fi ; }

while getopts ":c:d:a:p:n:h" opt; do
    case $opt in
        c) 
            structureFile=$OPTARG
            ;;
        d)
            mdp_dir=$OPTARG
            ;;
        p)  
            plumedFile=$OPTARG
            ;;
        h)
            usage
            ;;
        a)
            maxwarn=$OPTARG
            ;;
        n)
            indexFile=$OPTARG
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
if [ -z "$structureFile" ] || [ -z "$mdp_dir" ]; then
    usage
fi

# optional maxwarn not specified
if [ -z "$maxwarn" ]; then
    maxwarn=0
fi

structureName=`echo $structureFile | rev | cut -f1 -d"/" | rev | cut -f1 -d"."`

if [ -f nvt/"$structureName".xtc ] & [ -f npt/"$structureName".xtc ]; then
    read -p "An existing equilibration has been found. Do you want to use it for this run? (n)" -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        run_equilibration="true"
    else
        run_equilibration="false"
    fi
else
    run_equilibration="true"
fi


if [ "$run_equilibration" == "true" ]; then
    # temperature equilibration
    mkdir nvt
    gmx grompp -f "$mdp_dir"/nvt.mdp -c em/"$structureName".gro -r em/"$structureName".gro -p "$structureName".top -o nvt/"$structureName".tpr -po nvt/"$structureName".mdp -maxwarn "$maxwarn" -n "$indexFile" || { echo "Error: grompp for nvt failed" ; cleanup nvt ; exit 1; }
    gmx mdrun -v -s nvt/"$structureName".tpr -c nvt/"$structureName".gro -x nvt/"$structureName".xtc -cpo nvt/"$structureName.cpt" -e nvt/"$structureName".edr -g nvt/"$structureName".log || { echo "-> Error: gmx mdrun for nvt failed" ; cleanup nvt ; exit 1; }
    cleanup nvt

    # pressure equilibration
    mkdir npt
    gmx grompp -f "$mdp_dir"/npt.mdp -c nvt/"$structureName".gro -r nvt/"$structureName".gro -t nvt/"$structureName".cpt -p "$structureName".top -o npt/"$structureName".tpr -po npt/"$structureName".mdp -maxwarn "$maxwarn" -n "$indexFile" || { echo "Error: grompp for npt failed" ; cleanup npt ; exit 1; }
    gmx mdrun -v -s npt/"$structureName".tpr -c npt/"$structureName".gro -x npt/"$structureName".xtc -cpo npt/"$structureName.cpt" -e npt/"$structureName".edr -g npt/"$structureName".log || { echo "-> Error: gmx mdrun for npt failed" ; cleanup npt ; exit 1; }
    cleanup npt
fi


# production run
mkdir md0
gmx grompp -f "$mdp_dir"/md0.mdp -c npt/"$structureName".gro -t npt/"$structureName".cpt -p "$structureName".top -o md0/"$structureName".tpr -po md0/"$structureName".mdp -maxwarn "$maxwarn" -n "$indexFile" || { echo "Error: grompp for md0 failed" ; cleanup md0 ; exit 1; }

if [ ! -z "$plumedFile" ]; then
    gmx mdrun -v -s md0/"$structureName".tpr -c md0/"$structureName".gro -x md0/"$structureName".xtc -cpo md0/"$structureName".cpt -e md0/"$structureName".edr -g md0/"$structureName".log -plumed "$plumedFile" || { echo "-> Error: gmx mdrun for md0 failed" ; cleanup md0 ; exit 1; }
else
    gmx mdrun -v -s md0/"$structureName".tpr -c md0/"$structureName".gro -x md0/"$structureName".xtc -cpo md0/"$structureName".cpt -e md0/"$structureName".edr -g md0/"$structureName".log || { echo "-> Error: gmx mdrun for md0 failed" ; cleanup md0 ; exit 1; }
fi
cleanup md0
