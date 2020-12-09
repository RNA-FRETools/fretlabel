#!/bin/bash

usage() { echo "Run a multi replicate MD simulation
Usage: multi_run.sh -c <structure file (.gro)> -d <mdp directory> -a <max. warnings for grompp> (optional) -p <plumed file (.dat)> -n <index file (.ndx, optional)> -m <number of replicates (default: 8)> -e <run only this ensemble (nvt, np, md), default: run all>" 1>&2; exit 1; }
invalidOpt() { echo "Invalid option: -$OPTARG" 1>&2; exit 1; }
missingArg() { echo "Option -$OPTARG requires an argument" 1>&2; exit 1; }
cleanup() { if ls -f $1/\#* 1> /dev/null 2>&1 ; then rm $1/\#* ; fi ; }

while getopts ":c:d:a:p:n:m:e:h" opt; do
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
        m)
            replicates=$OPTARG
            ;;
        e)
            ensemble=$OPTARG
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

# optional ensemble not specified
if [ -z "$ensemble" ]; then
    ensemble="all"
fi

structureName=`echo $structureFile | rev | cut -f1 -d"/" | rev | cut -f1 -d"."`
folders=`seq -s' ' 0 $replicates`

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
    if [ "$ensemble" == "all" ] || [ "$ensemble" == "nvt" ] ; then
    mkdir nvt
    for i in $folders; do
        gmx grompp -f "$mdp_dir"/nvt.mdp -c em/"$structureName".gro -r em/"$structureName".gro -p "$structureName".top -o nvt/$i/"$structureName".tpr -po nvt/$i/"$structureName".mdp -maxwarn "$maxwarn" -n "$indexFile" || { echo "Error: grompp for nvt failed" ; cleanup nvt ; exit 1; }
    done
    mpirun -np $replicates gmx mdrun -v -s nvt/$i/"$structureName".tpr -c nvt/$i/"$structureName".gro -x nvt/$i/"$structureName".xtc -cpo nvt/$i/"$structureName.cpt" -e nvt/$i/"$structureName".edr -g nvt/$i/"$structureName".log -multidir $folders || { echo "-> Error: gmx mdrun for nvt failed" ; cleanup nvt ; exit 1; }
    cleanup nvt
    fi

    # pressure equilibration
    if [ "$ensemble" == "all" ] || [ "$ensemble" == "npt" ] ; then
    mkdir npt
    for i in {0..$replicates}; do
        gmx grompp -f "$mdp_dir"/npt.mdp -c nvt/$i/"$structureName".gro -r nvt/$i/"$structureName".gro -t nvt/$i/"$structureName".cpt -p "$structureName".top -o npt/$i/"$structureName".tpr -po npt/$i/"$structureName".mdp -maxwarn "$maxwarn" -n "$indexFile" || { echo "Error: grompp for npt failed" ; cleanup npt ; exit 1; }
    done
    mpirun -np $replicates gmx mdrun -v -s npt/$i/"$structureName".tpr -c npt/$i/"$structureName".gro -x npt/$i/"$structureName".xtc -cpo npt/$i/"$structureName.cpt" -e npt/$i/"$structureName".edr -g npt/$i/"$structureName".log -multidir $folders || { echo "-> Error: gmx mdrun for npt failed" ; cleanup npt ; exit 1; }
    cleanup npt
    fi
fi


# production run
mkdir md0
if [ "$ensemble" == "all" ] || [ "$ensemble" == "md" ] ; then
for i in {0..$replicates}; do
    gmx grompp -f "$mdp_dir"/md0.mdp -c npt/$i/"$structureName".gro -t npt/$i/"$structureName".cpt -p "$structureName".top -o md0/$i/"$structureName".tpr -po md0/$i/"$structureName".mdp -maxwarn "$maxwarn" -n "$indexFile" || { echo "Error: grompp for md0 failed" ; cleanup md0 ; exit 1; }
done
if [ ! -z "$plumedFile" ]; then
    mpirun -np $replicates gmx mdrun -v -s md0/$i/"$structureName".tpr -c md0/$i/"$structureName".gro -x md0/$i/"$structureName".xtc -cpo md0/$i/"$structureName".cpt -e md0/$i/"$structureName".edr -g md0/$i/"$structureName".log -plumed "$plumedFile" -multidir $folders || { echo "-> Error: gmx mdrun for md0 failed" ; cleanup md0 ; exit 1; }
else
    mpirun -np $replicates gmx mdrun -v -s md0/$i/"$structureName".tpr -c md0/$i/"$structureName".gro -x md0/$i/"$structureName".xtc -cpo md0/$i/"$structureName".cpt -e md0/$i/"$structureName".edr -g md0/$i/"$structureName".log -multidir $folders || { echo "-> Error: gmx mdrun for md0 failed" ; cleanup md0 ; exit 1; }
fi
cleanup md0
fi