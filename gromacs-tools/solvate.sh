#!/bin/bash

# Solvation of an biomolecule in a water box for MD simulations

# cmd parsing functions
usage() { echo "Solvate of an biomolecule in a water box for an MD simulation
Usage: solvate.sh -f <structure file>" 1>&2; exit 1; }
invalidOpt() { echo "Invalid option: -$OPTARG" 1>&2; exit 1; }
missingArg() { echo "Option -$OPTARG requires an argument" 1>&2; exit 1; }


#------------
# cmd parsing
#------------

while getopts ":s:h" opt; do
    case $opt in
        f) 
            structureFile=$OPTARG
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
if [ -z "$structureFile" ]; then
    usage
fi

# check if file is present
if [ ! -e "$structureFile" ]; then
    echo Error: the specified stucture file does not exist
    exit 1
fi

# GROMACS pipeline

structureName=`echo $structureFile | cut -f1 -d"."`
gmx pdb2gmx -f "$structureFile" -o "$structureName".gro -p "$structureName".top -i "$structureName"_posre.itp -merge all -ter yes || exit 1
gmx editconf -f "$structureName".gro -o "$structureName"_box.gro -bt dodecahedron -d 1 || exit 1
gmx solvate -cp "$structureName"_box.gro -o "$structureName"_solv.gro -p "$structureName".top || exit 1
gmx grompp -f ions.mdp -c "$structureName"_solv.gro -p "$structureName".top -o "$structureName".tpr -po "$structureName".mdp || exit 1
gmx genion -s "$structureName".tpr -o "$structureName"_withIons.gro -p "$structureName".top -nname Cl -pname K -neutral || exit 1
gmx grompp -f em.mdp -c "$structureName"_withIons.gro -p "$structureName".top -o "$structureName"_min.tpr -po "$structureName"_min.mdp || exit 1
gmx mdrun -s "$structureName"_min || exit 1

# clean up
if ls -f \#* 1> /dev/null 2>&1 ; then
    rm \#*
fi
