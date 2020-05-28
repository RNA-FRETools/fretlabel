#!/bin/bash

# Solvation of an biomolecule in a water box for MD simulations

# cmd parsing functions
usage() { echo "Run a single replicate MD simulation
Usage: single_run.sh -c <structure file (.gro)> -d <mdp directory> -a <max. warnings for grompp> (optional)" 1>&2; exit 1; }
invalidOpt() { echo "Invalid option: -$OPTARG" 1>&2; exit 1; }
missingArg() { echo "Option -$OPTARG requires an argument" 1>&2; exit 1; }
cleanup() { if ls -f $1/\#* 1> /dev/null 2>&1 ; then rm $1/\#* ; fi ; }

#------------
# cmd parsing
#------------

while getopts ":c:d:a:h" opt; do
    case $opt in
        c) 
            structureFile=$OPTARG
            ;;
        d)
            mdp_dir=$OPTARG
            ;;
        h)
            usage
            ;;
        a)
            maxwarn=$OPTARG
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

structureName=`echo $structureFile | rev | cut -f1 -d"/" | rev | cut -f1 -d"."`

# temperature equilibration
mkdir nvt
gmx grompp -f "$mdp_dir"/nvt.mdp -c em/"$structureName".gro -r em/"$structureName".gro -p "$structureName".top -o nvt/"$structureName".tpr -po nvt/"$structureName".mdp -maxwarn "$maxwarn" || { echo "Error: grompp for nvt failed" ; cleanup nvt ; exit 1; }

gmx mdrun -v -s nvt/"$structureName".tpr -c nvt/"$structureName".gro -x nvt/"$structureName".xtc -cpo nvt/"$structureName.cpt" -e nvt/"$structureName".edr -g nvt/"$structureName".log || { echo "-> Error: gmx mdrun for nvt failed" ; cleanup nvt ; exit 1; }

cleanup nvt


# pressure equilibration
mkdir npt
gmx grompp -f "$mdp_dir"/npt.mdp -c nvt/"$structureName".gro -r nvt/"$structureName".gro -t nvt/"$structureName".cpt -p "$structureName".top -o npt/"$structureName".tpr -po npt/"$structureName".mdp -maxwarn "$maxwarn" || { echo "Error: grompp for npt failed" ; cleanup npt ; exit 1; }

gmx mdrun -v -s npt/"$structureName".tpr -c npt/"$structureName".gro -x npt/"$structureName".xtc -cpo npt/"$structureName.cpt" -e npt/"$structureName".edr -g npt/"$structureName".log || { echo "-> Error: gmx mdrun for npt failed" ; cleanup npt ; exit 1; }

cleanup npt


# production run
mkdir md0
gmx grompp -f "$mdp_dir"/md0.mdp -c npt/"$structureName".gro -t npt/"$structureName".cpt -p "$structureName".top -o md0/"$structureName".tpr -po md0/"$structureName".mdp -maxwarn "$maxwarn" || { echo "Error: grompp for md0 failed" ; cleanup md0 ; exit 1; }

gmx mdrun -v -s md0/"$structureName".tpr -c md0/"$structureName".gro -x md0/"$structureName".xtc -cpo md0/"$structureName".cpt -e md0/"$structureName".edr -g md0/"$structureName".log || { echo "-> Error: gmx mdrun for md0 failed" ; cleanup md0 ; exit 1; }

cleanup md0
