#!/bin/bash

# Solvation of an biomolecule in a water box for MD simulations

# cmd parsing functions
usage() { echo "Run a single replicate MD simulation
Usage: single_run.sh -f <structure file (.gro, .pdb)> -d <mdp directory>" 1>&2; exit 1; }
invalidOpt() { echo "Invalid option: -$OPTARG" 1>&2; exit 1; }
missingArg() { echo "Option -$OPTARG requires an argument" 1>&2; exit 1; }
cleanup() { if ls -f $1/\#* 1> /dev/null 2>&1 ; then rm $1/\#* ; fi ; exit 1; }

#------------
# cmd parsing
#------------

while getopts ":c:d:h" opt; do
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

mkdir nvt
gmx grompp -f "$mdp_dir"/nvt.mdp -c em/"$structureFile""$i".gro -p nvt/"$structureFile".top -o nvt/"$structureFile".tpr -po nvt/"$structureFile".mdp || { echo "Error: grompp for nvt failed" ; cleanup "nvt"; }

mkdir npt
gmx grompp -f "$mdp_dir"/nvt.mdp -c nvt/"$structureFile""$i".gro -t nvt/"$inputName""$i".cpt -p npt/"$structureFile".top -o npt/"$structureFile".tpr -po npt/"$structureFile".mdp || { echo "Error: grompp for nvt failed" ; cleanup "npt"; }