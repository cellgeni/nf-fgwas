#!/bin/bash

#######################################################################
#                  submit job for prepare_sum_stat                    #
#######################################################################

while IFS="" read -r var || [ -n "$var" ]
do
	if [ -z "$var" ]; then continue; fi
	if [[ ${var:0:1} == "#" ]]
	then
		echo "skip line: $var"
	else
		study=$(echo "$var" | cut -d',' -f1)
		file_path=$(echo "$var" | cut -d',' -f2)
		echo "file: $file_path"

		runID=${study}_$(date +"%d%b%y")

		# prepare output folder
		OUTDIR="prep_sum_stat/${runID}"

		mkdir -p "$OUTDIR/FarmOut"
		rm -f \
			"$OUTDIR/FarmOut/stderr.log" \
			"$OUTDIR/FarmOut/stdout.log"

		# submit job
		echo "submitting job for $runID..."

		MEM=50G
		bsub \
			-q small \
			-J"prepare_inputsJob_$runID" \
			-e "$OUTDIR/FarmOut/stderr.log" \
			-o "$OUTDIR/FarmOut/stdout.log" \
			-n 8 \
			-M"$MEM" \
			-R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" \
			"singularity exec --bind /lustre:/lustre /nfs/cellgeni/singularity/images/fgwas_v0-2.sif python nf-gwas/bin/check_summ_stat.py -i $file_path -o $OUTDIR 1>&2"

	fi
done < <(grep -v "^$" summ_stat_paths.txt)

# # watch
watch -n 60 "bjobs -w | grep prepare_inputs | grep PEND | wc -l | xargs printf 'number of pending jobs: %s\n'; bjobs -w | grep prepare_inputs | grep RUN | wc -l | xargs printf 'number of running jobs: %s\n\nPress Ctrl-C to stop watching.\n'"

# alternative for single job
# echo "PEND ...waiting for job to start"
# while [ -n "$(bjobs -w | grep $runID | grep PEND)" ]
# do
# 	sleep 10s
# done
# if [ -n "$(bjobs -w | grep $runID | grep RUN)" ]
# then
# 	echo "RUN ...showing stderr:"
# 	touch "$OUTDIR/FarmOut/stderr.log"
# 	tail -f "$OUTDIR/FarmOut/stderr.log"
# else
# 	echo "terminated ...showing stdout:"
# 	tail -f "$OUTDIR/FarmOut/stdout.log"
# fi
