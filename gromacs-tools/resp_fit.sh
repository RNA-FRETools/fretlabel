#!/bin/bash

usage() { echo "Two-stage restrained electrostatic potential (RESP) fitting
Usage: resp_fit.sh -n <name if the input mol2 file (no extension)> -i <input directory> -o <output directory> -g <capping group file (.dat)> -c <net charge>" 1>&2; exit 1; }
invalidOpt() { echo "Invalid option: -$OPTARG" 1>&2; exit 1; }
missingArg() { echo "Option -$OPTARG requires an argument" 1>&2; exit 1; }


while getopts ":n:i:o:g:c:h" opt; do
    case $opt in
        n) 
            name=$OPTARG
            ;;
        i)
            input_folder=$OPTARG
            ;;
        o)  
            output_folder=$OPTARG
            ;;
        g)
            capping_group=$OPTARG
            ;;
        c)
            net_charge=$OPTARG
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

# any cmd line argument not given
if [ -z "$name" ] || [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$capping_group" ] || [ -z "$net_charge" ]; then
    usage
fi

antechamber -i "$input_folder/$name".mol2 -fi mol2 -o "$output_folder/$name".ac -fo ac -pf yes -nc "$net_charge" || exit 1

n_atom=`awk '$1 == "GROUP" {print $2}' "$input_folder/$capping_group"`
group_constraint=`awk '$1 == "GROUP" {print $3}' "$input_folder/$capping_group"`

respgen -i "$output_folder/$name".ac -o "$output_folder/$name".respin1 -f resp1 -a "$input_folder/$capping_group" || exit 1
respgen -i "$output_folder/$name".ac -o "$output_folder/$name".respin2 -f resp2 -a "$input_folder/$capping_group" || exit 1

# since respgen rounds the group constraint to three decimals replace it with the value from the capping group
for i in `seq $(echo $n_atom | wc -w)`;do
    atom=`echo $n_atom | cut -f$i -d' '`
    group=`echo $group_constraint | cut -f$i -d' '`
    sed -z -i "s/$atom\s*-*\w\.\w*$p/$atom  $group/$i" "$output_folder/$name".respin1 # -i (inplace) -z (treat file as single line)
    sed -z -i "s/$atom\s*-*\w\.\w*$p/$atom  $group/$i" "$output_folder/$name".respin2
done || exit 1

espgen -i "$output_folder/$name"_hf_esp.gout -o "$output_folder/$name"_hf_esp.esp || exit 1
mv QIN "$output_folder"/

resp -O -i "$output_folder/$name".respin1 -o "$output_folder/$name".respout1 -e "$output_folder/$name"_hf_esp.esp -q  "$output_folder"/QIN -t "$output_folder"/qout_stage1 -p "$output_folder"/punch1 -s "$output_folder"/esout1 || exit 1
resp -O -i "$output_folder/$name".respin2 -o "$output_folder/$name".respout2 -e "$output_folder/$name"_hf_esp.esp -q "$output_folder"/qout_stage1 -t "$output_folder"/qout_stage2 -p "$output_folder"/punch2 -s "$output_folder"/esout2 || exit 1

antechamber -i "$output_folder/$name".ac -fi ac -o "$output_folder/$name"_resp.mol2 -fo mol2 -c rc -cf "$output_folder"/qout_stage2 -pf yes -at amber
