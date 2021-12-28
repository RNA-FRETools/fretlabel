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
folders=`seq -s' ' 1 $replicates`

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
    nvt_folders=`for i in $(seq -s' ' 1 $replicates); do echo npt/$i;done`
    for i in $folders; do
        mkdir $i
        gmx grompp -f "$mdp_dir"/nvt.mdp -c em/"$structureName".gro -r em/"$structureName".gro -p "$structureName".top -o nvt/$i/"$structureName".tpr -po nvt/$i/"$structureName".mdp -maxwarn "$maxwarn" -n "$indexFile" || { echo "Error: grompp for nvt failed" ; cleanup nvt/$i ; exit 1; }
    done
    mpirun -np $replicates gmx mdrun -v -s "$structureName".tpr -c "$structureName".gro -x "$structureName".xtc -cpo "$structureName.cpt" -e "$structureName".edr -g "$structureName".log -multidir $nvt_folders || { echo "-> Error: gmx mdrun for nvt failed" ; for i in $folders; do cleanup nvt/$i;done ; exit 1; }
    for i in $folders; do
        cleanup nvt/$i
    done
    fi

    # pressure equilibration
    if [ "$ensemble" == "all" ] || [ "$ensemble" == "npt" ] ; then
    mkdir npt
    npt_folders=`for i in $(seq -s' ' 1 $replicates); do echo npt/$i;done`
    for i in $folders; do
        mkdir npt/$i
        gmx grompp -f "$mdp_dir"/npt.mdp -c nvt/$i/"$structureName".gro -r nvt/$i/"$structureName".gro -t nvt/$i/"$structureName".cpt -p "$structureName".top -o npt/$i/"$structureName".tpr -po npt/$i/"$structureName".mdp -maxwarn "$maxwarn" -n "$indexFile" || { echo "Error: grompp for npt failed" ; cleanup npt/$i ; exit 1; }
    done
    mpirun -np $replicates gmx mdrun -v -s "$structureName".tpr -c "$structureName".gro -x "$structureName".xtc -cpo "$structureName.cpt" -e "$structureName".edr -g "$structureName".log -multidir $npt_folders || { echo "-> Error: gmx mdrun for npt failed" ; for i in $folders; do cleanup npt/$i;done ; exit 1; }
    for i in $folders; do
        cleanup npt/$i
    done
    fi
fi


# production run
if [ "$ensemble" == "all" ] || [ "$ensemble" == "md" ] ; then
mkdir md0
for i in $folders; do
    mkdir md0/$i
    md0_folders=`for i in $(seq -s' ' 1 $replicates); do echo md0/$i;done`
    gmx grompp -f "$mdp_dir"/md0.mdp -c npt/$i/"$structureName".gro -t npt/$i/"$structureName".cpt -p "$structureName".top -o md0/$i/"$structureName".tpr -po md0/$i/"$structureName".mdp -maxwarn "$maxwarn" -n "$indexFile" || { echo "Error: grompp for md0 failed" ; cleanup md0/$i ; exit 1; }
done
if [ ! -z "$plumedFile" ]; then
    mpirun -np $replicates gmx mdrun -v -s "$structureName".tpr -c "$structureName".gro -x "$structureName".xtc -cpo "$structureName".cpt -e "$structureName".edr -g "$structureName".log -plumed "$plumedFile" -multidir $md0_folders || { echo "-> Error: gmx mdrun for md0 failed" ; for i in $folders; do cleanup md0/$i;done ; exit 1; }
else
    mpirun -np $replicates gmx mdrun -v -s "$structureName".tpr -c "$structureName".gro -x "$structureName".xtc -cpo "$structureName".cpt -e "$structureName".edr -g "$structureName".log -multidir $md0_folders || { echo "-> Error: gmx mdrun for md0 failed" ; for i in $folders; do cleanup md0/$i;done ; exit 1; }
fi
for i in $folders; do
    cleanup md0/$i
done
fi
