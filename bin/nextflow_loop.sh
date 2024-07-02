#!/bin/bash

set -eo pipefail

# module load cellgen/singularity
# module load cellgen/nextflow

LDIR="$(pwd)"

cell_type_file=$(readlink -f "29Jun24_prepare_inputs/tss_cell_type_exp.txt")
tss_file="${cell_type_file}.gz"

while IFS="" read -r var || [ -n "$var" ]
do
	if [ -z "$var" ]; then continue; fi
	if [[ ${var:0:1} == "#" ]]
	then
		echo "skip line: $var"
	else
		study=$(echo "$var" | cut -d',' -f1)
		gwas_path=$(echo "$var" | cut -d',' -f2)

		echo ">RUN line: ID: $study file: $gwas_path"

		runID=${study}_$(date +"%d%b%y")

		# prepare output folder
		OUTDIR="$LDIR/fGWAS_results/${runID}"

		screen -S "$runID" -dm bash -c "/bin/bash $LDIR/scripts/run_nextflow.sh $OUTDIR $LDIR $runID $tss_file $cell_type_file $gwas_path"
	fi
done < <(grep -v "^$" summ_stat_paths_fixed.txt)

