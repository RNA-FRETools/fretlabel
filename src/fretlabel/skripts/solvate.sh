#!/bin/bash

# cmd parsing functions
usage() { echo "Solvate of an biomolecule in a water box for an MD simulation
Usage: solvate.sh -f <structure file (.gro, .pdb)> -w <solvent file (tip4p.gro, spc216.gro)> -d <mdp directory> -b <distance to box edge (default: 1nm)>" 1>&2; exit 1; }
invalidOpt() { echo "Invalid option: -$OPTARG" 1>&2; exit 1; }
missingArg() { echo "Option -$OPTARG requires an argument" 1>&2; exit 1; }
cleanup() { if ls -f $1/\#* 1> /dev/null 2>&1 ; then rm $1/\#* ; fi ; }


while getopts ":f:w:d:b:h" opt; do
    case $opt in
        f) 
            structureFile=$OPTARG
            ;;
        d)
            mdp_dir=$OPTARG
            ;;
        w)  
            waterFile=$OPTARG
            ;;
        b) 
            dist2boxedge=$OPTARG
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
if [ -z "$structureFile" ] || [ -z "$mdp_dir" ] || [ -z "$waterFile" ]; then
    usage
fi

# check if file is present
if [ ! -e "$structureFile" ]; then
    echo Error: the specified structure file does not exist
    exit 1
fi

if [ -z $dist2boxedge ]; then
    dist2boxedge=1.5
fi

# note: .top must be in root directory
structureName=`echo $structureFile | cut -f1 -d"."`
mkdir em

gmx pdb2gmx -f "$structureFile" -o em/"$structureName".gro -p "$structureName".top -i em/"$structureName".itp || { echo "-> Error: gmx pdb2gmx failed" ; cleanup em; exit 1; }

gmx editconf -f em/"$structureName".gro -o em/"$structureName".gro -bt dodecahedron -d $dist2boxedge || { echo "-> Error: gmx editconf failed" ; cleanup em; exit 1; }

gmx solvate -cp em/"$structureName".gro -cs "$waterFile" -o em/"$structureName".gro -p "$structureName".top || { echo "Error: gmx solvate failed" ; cleanup em; exit 1; }

gmx grompp -f "$mdp_dir"/em.mdp -c em/"$structureName".gro -p "$structureName".top -o em/"$structureName".tpr -po em/"$structureName".mdp -maxwarn 2 || { echo "-> Error: grompp after solvation failed" ; cleanup em; exit 1; }

echo SOL | gmx genion -s em/"$structureName".tpr -o em/"$structureName".gro -p "$structureName".top -nname Cl -pname K -neutral || { echo "-> Error: gmx genion failed" ; cleanup em; exit 1; }

gmx grompp -f "$mdp_dir"/em.mdp -c em/"$structureName".gro -p "$structureName".top -o em/"$structureName".tpr -po em/"$structureName".mdp || { echo "-> Error: gmx grompp after genion failed" ; cleanup em; exit 1; }

gmx mdrun -v -s em/"$structureName".tpr -c em/"$structureName".gro -o em/"$structureName".trr -e em/"$structureName".edr -g em/"$structureName".log || { echo "-> Error: gmx mdrun failed" ; cleanup em; exit 1; }

gmx editconf -f em/"$structureName".gro -o em/"$structureName".pdb || { echo "-> Error: gmx editconf failed" ; cleanup em; exit 1; }

cleanup ./
cleanup em
