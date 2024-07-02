#!/bin/bash

#######################################################################
#                  submit job for prepare_inputs.py                   #
#######################################################################

# ID of run instance (e.g. loop over to submit several jobs)
runID=$(date +"%d%b%y")

# prepare output folder
OUTDIR="${runID}_prepare_inputs"

mkdir -p "$OUTDIR/FarmOut"
rm -f \
	"$OUTDIR/FarmOut/stderr.log" \
	"$OUTDIR/FarmOut/stdout.log"

# submit job
echo "submitting job for $runID..."

MEM=50G
bsub \
	-q small \
	-J"prepare_inputsJob$runID" \
	-e "$OUTDIR/FarmOut/stderr.log" \
	-o "$OUTDIR/FarmOut/stdout.log" \
	-n 8 \
	-M"$MEM" \
	-R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" \
	"singularity exec --bind /lustre:/lustre /nfs/cellgeni/singularity/images/fgwas_v0-2.sif python nf-gwas/bin/prepare_input.py tss -i $1 -g $2 -o $OUTDIR 1>&2"

# # watch
# watch -n 60 "bjobs -w | grep prepare_inputs | grep PEND | wc -l | xargs printf 'number of pending jobs: %s\n'; bjobs -w | grep prepare_inputs | grep RUN | wc -l | xargs printf 'number of running jobs: %s\n\nPress Ctrl-C to stop watching.\n'"
# alternative for single job
echo "PEND ...waiting for job to start"
while [ -n "$(bjobs -w | grep $runID | grep PEND)" ]
do
	sleep 10s
done
if [ -n "$(bjobs -w | grep $runID | grep RUN)" ]
then
	echo "RUN ...showing stderr:"
	tail -f "$OUTDIR/FarmOut/stderr.log"
else
	echo "terminated ...showing stdout:"
	tail -f "$OUTDIR/FarmOut/stdout.log"
fi
