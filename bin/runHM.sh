#!/bin/bash

#######################################################################
#                          HANDLE PARAMETERS                          #
#######################################################################

# input ldsc_input: <getLDSC/collected/input_hm.gz>
# input rna_infile: <tmp/tss_cell_type_exp.txt.gz>
# input broad_fine_mapping: <input/broad_fine_mapping.tsv>
# output hm_done: <done/hm.done>
# output outdir: <HM>

LDSC_INPUT=""
ATAC_CELL_TYPES="tmp/ATAC_peaks_per_celltype_wide_sorted_merged.tsv"
RNA_CELL_TYPES="input/tss_cell_type_exp.txt"
RNA_INFILE="tmp/tss_cell_type_exp.txt.gz"
BROAD_FINE_MAPPING="input/broad_fine_mapping.tsv"

PARAMS=""   # positional arguments: STUDY

ATAC=""     # include atac data flag
JOBIND=""   # jobindex (corresponds to cell type)
OUTDIR="HM" # name of output directory
NROW=""     # number of rows
NCELLT=""   # number of cell types

while (( "$#" ))
do
	case "$1" in
		-h|--help)
			# show help
			shift
			;;
		-a|--atac)
			# covid flag
			ATAC="--atac"
			shift
			;;
		-j|--jobindex)
			# set jobindex
			if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
				JOBIND="$2"
				shift 2
			else
				echo "Error: Argument for $1 is missing" >&2
				exit 1
			fi
			;;
		-i|--ldscinput)
			# set output directory
			if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
				LDSC_INPUT="$2"
				shift 2
			else
				echo "Error: Argument for $1 is missing" >&2
				exit 1
			fi
			;;
		-o|--outdir)
			# set output directory
			if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
				OUTDIR="$2"
				shift 2
			else
				echo "Error: Argument for $1 is missing" >&2
				exit 1
			fi
			;;
		-n|--nrow)
			# set number of rows
			if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
				NROW="$2"
				shift 2
			else
				echo "Error: Argument for $1 is missing" >&2
				exit 1
			fi
			;;
		-c|--ncellt)
			# set number of cell types
			if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
				NCELLT="$2"
				shift 2
			else
				echo "Error: Argument for $1 is missing" >&2
				exit 1
			fi
			;;
		-*|--*=) 
			# unsupported flags
			echo "Error: Unsupported flag $1" >&2
			exit 1
			;;
		*) 
			# positional arguments
			PARAMS="$PARAMS $1"
			shift
			;;
	esac
done
PARGS=($PARAMS)

STUDY=${PARGS[0]}

# check required params / show help
if [ -z "$STUDY" ]
then
	printf -- "\nusage: runHM.sh STUDY [--atac] [--nrow NR] [--ncellt NC]\n\n"
	printf -- "STUDY        ID of study\n\n"
	printf -- "-a, --atac   include ATAC data in analysis\n\n"
	printf -- "the following arguments will be set automatically:\n"
	printf -- "-j, -jobindex JI  number of rows in input_hm.gz file\n"
	printf -- "-i, -ldscinput PATH   path to the input_hm.gz file\n"
	printf -- "-o, -outdir DIR   name of output directory (default: HM)\n"
	printf -- "-n, --nrow NR     number of rows in input_hm.gz file\n"
	printf -- "-c, --ncellt NC   number of cell types in $RNA_INFILE\n\n"
	exit
fi

if [ -z "$NROW" ]
then
	NROW=$(zcat "$LDSC_INPUT" | wc | awk '{print $1}')
fi

if [ -z "$NCELLT" ]
then
	NCELLT=`zcat $RNA_INFILE | head -n 1 | wc | awk '{print $2-2}'`
fi

#######################################################################
#                             submit job                              #
#######################################################################


if [ -z "$JOBIND" ]
then
	#NROW=`bash /software/team170/nk5/fGWAS/collect.sh $STUDY WOMHC`
	#NROW=`bash /software/team170/nk5/fGWAS/collect.sh $STUDY`

	mkdir -p $STUDY/HM
	mkdir -p $STUDY/FarmOut_runHM

	rm -f -r $STUDY/HM/*
	rm -f $STUDY/FarmOut_runHM/*

	MEM=20000

	temp_file=$(mktemp)
	trap "rm -f $temp_file" 0 2 3 15

	# submit job
	#-w "done(LDSCJob$STUDY)" \
	bsub -q long -m modern_hardware -J"hmJob[1-$NCELLT]" -e $STUDY/FarmOut_runHM/stderr.%I -o $STUDY/FarmOut_runHM/stdout.%I -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" "bash runHM.sh $STUDY $ATAC --jobindex \$LSB_JOBINDEX --nrow $NROW --ncellt $NCELLT" | tee $temp_file
	
	jobid=$(grep -Po "(?<=Job <)(\d+)(?=>)" $temp_file)

	# watch job
	watch -n 60 "bjobs -w | grep hmJob | grep $jobid | grep PEND | wc -l | xargs printf '"$(pwd)"\n\n> jobs for hmJob $STUDY $ATAC\n\nnumber of pending jobs: %s\n'; bjobs -w | grep hmJob | grep $jobid | grep RUN | wc -l | xargs printf 'number of running jobs: %s\n\nPress Ctrl-C to stop watching.\n'"

	exit
fi

#######################################################################
#                               run job                               #
#######################################################################

printf -- "\n--------------------\n"
printf -- "STUDY: $STUDY\n"
printf -- "ATAC: $ATAC\n"
printf -- "NROW: $NROW\n"
printf -- "NCELLT: $NCELLT\n"
printf -- "JOBIND: $JOBIND\n"
printf -- "all params: $PARAMS\n"
printf -- "--------------------\n\n"

##############################
#  column line for TSS file  #
##############################

NFE=`zcat $RNA_INFILE | wc | awk '{print $1}'`
NCOL=`zcat $RNA_INFILE | head -n 1 | wc | awk '{print $2}'`
echo $NCOL
COL=`echo $JOBIND | awk -v TOT=$NCOL '
	{
		printf "S,S"; 
		for(i=3;i<TOT;i++){
			if(i-2==$1){
				printf ",N0";
			}else{
				printf ",S"
			}
		}; 
		printf ",N0";}'`

###################################
#  column line for input_hm file  #
###################################

NCOL2=`zcat "$LDSC_INPUT" | head -n 1 | wc | awk '{print $2}'`

if [ -z "$ATAC" ]
then
	NUM_COLS=`expr $NCOL2 '-' 3`
	[ $NUM_COLS -eq 0 ] && ADD_COLS="" || ADD_COLS=$(printf ',S%.0s' $(seq 1 $NUM_COLS))
	COL2="I,B,O"$ADD_COLS
else
	NCOL=`zcat $RNA_INFILE | head -n 1 | wc | awk '{print $2}'`
	# get rough cell type category for ATAC instead of TSS file cell type ($JOBIND)
	# the mapping between fine and rough cell type annotation is read from
	# the BED file `ATAC_peaks_per_celltype.bed` and translated forth and back
	# from/to column numbers of the other input files (tss and input_hm)
	TSS_FILE_CELL_TYPES=`cat $RNA_CELL_TYPES | head -n 1 | cut -d$'\t' --complement -f"1,2,$NCOL"`
	ATAC_FILE_CELL_TYPES=`cat $ATAC_CELL_TYPES | head -n 1 | cut -d$'\t' --complement -f1`

	COL2=`cat $BROAD_FINE_MAPPING | awk -v TOT="$NCOL2" -v JOBIND="$JOBIND" -v TSS_CT="$TSS_FILE_CELL_TYPES" -v ATAC_CT="$ATAC_FILE_CELL_TYPES" '
		BEGIN{
			# mappings from/to column numbers
			split(TSS_CT, num2cellt_tss, "\t");
			split(ATAC_CT, num2cellt_atac, "\t");
			for(i=1; i<=length(num2cellt_atac); i++){
				cellt = num2cellt_atac[i];
				cellt2num_atac[cellt] = i;
			}
		}
		{
			# mapping between fine and rough cell type categories
			cellt_map[$1] = $2;
		}
		END{
			# write column line
			printf "I,B,O"; 
			job_ct = num2cellt_tss[JOBIND];
			atac_ct = cellt_map[job_ct];
			atac_ct_num = cellt2num_atac[atac_ct];
			for(i=4;i<=TOT;i++){
				if(i-3==atac_ct_num){
					printf ",C2";
					#printf ",S";
				}else{
					printf ",S";
				}
			};}'`
fi

#####################
#  execute command  #
#####################

# run_command="/nfs/users/nfs_n/nk5/Project/C/PHM/src/hm -v --ab-feature-level 1.0 100.0 -r $NROW -f $NFE -o $STUDY/$OUTDIR/Out$JOBIND -i $LDSC_INPUT -c $COL2 -j $RNA_INFILE -d $COL"

run_command="/nfs/team205/jp30/projects/code/PHM/src/hm -v --ab-feature-level 1.0 100.0 -r $NROW -f $NFE -o $STUDY/$OUTDIR/Out$JOBIND -i $LDSC_INPUT -c $COL2 -j $RNA_INFILE -d $COL"

echo "$run_command"

$run_command

