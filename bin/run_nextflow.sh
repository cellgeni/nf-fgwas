#!/bin/bash

. ~/.bashrc

set -x

module load cellgen/singularity
module load cellgen/nextflow

NX_FILE="/lustre/scratch127/cellgen/cellgeni/tickets/tic-3145/nf-gwas/main.nf"


OUTDIR="$1"
LDIR="$2"

runID="$3"
tss_file="$4"
cell_type_file="$5"
gwas_path="$6"


mkdir -p "$OUTDIR"
cd "$OUTDIR"

nextflow -C $LDIR/nextflow.config run $NX_FILE \
	--project_tag $runID \
	--tss_file $tss_file \
	--cell_types $cell_type_file \
	--gwas_path $gwas_path \

exec bash

