#!/bin/bash

GEX_INFILE="tss_cell_type_exp.txt.gz"
VCF_FILES_DIR="/nfs/cellgeni/pipeline-files/1000G/20190312/European/"

# The number of genes per job
NGENE=100
WSIZE=500000
PARQUET_FILE=""

PARAMS=""
IS_COVID=""
CUSTOM_GWAS=""
JOBIND=""
LIST=""


#######################################################################
#                          HANDLE PARAMETERS                          #
#######################################################################

while (( "$#" ))
do
	case "$1" in
		-h|--help)
			# show help
			shift
			;;
		-c|--covid)
			# covid flag
			IS_COVID="--covid"
			shift
			;;
		-l|--list)
			# list available studies
			LIST=true
			shift
			;;
		-g|--gwas)
			# set a custom GWAS study
			if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
				CUSTOM_GWAS="$2"
				shift 2
			else
				echo "Error: Argument for $1 is missing" >&2
				exit 1
			fi
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
		-n|--ngene)
			# set jobindex
			if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
				NGENE="$2"
				shift 2
			else
				echo "Error: Argument for $1 is missing" >&2
				exit 1
			fi
			;;
		-w|--windowsize)
			# set window size
			if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
				WSIZE="$2"
				shift 2
			else
				echo "Error: Argument for $1 is missing" >&2
				exit 1
			fi
			;;
		-p|--parquetfile)
			# set dir with parquet files
			if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
				PARQUET_FILE="$2"
				shift 2
			else
				echo "Error: Argument for $1 is missing" >&2
				exit 1
			fi
			;;
		-v|--vcffilesdir)
			# set dir with vcf files
			if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
				VCF_FILES_DIR="$2"
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
#eval set -- "$PARAMS"
PARGS=($PARAMS)
GWAS=${PARGS[0]}
ATAC=${PARGS[1]}

# list available studies
if [ -n "$LIST" ]
then
	ls /warehouse/cellgeni/otar-gwas-ss/gwas/ | while read f; do basename $f .parquet; done
	ls /nfs/users/nfs_n/nk5/sanger/GWAS/COVID19/Orig/ | while read f; do basename $f .txt.gz; done
	exit
fi


# check required params / show help
if [ -z "$GWAS" ]
then
	printf -- "\nusage: getLDSC.sh STUDY [ATAC]\n\n"
	printf -- "STUDY             ID of study\n"
	printf -- "ATAC              tabix indexed bed file with values of 0/1\n"
	printf -- "                  for each celltype per peak region\n\n"
	printf -- "-l, --list        list available study IDs and exit\n"
	printf -- "-g, --gwas        path to tsv file with summary statistics\n\n"
	printf -- "-n, --ngene       number of genes per chunk/batch\n\n"
	printf -- "-w, --windowsize  distance up and downstream of TSS to include\n"
	printf -- "                  SNPs (default 500k, i.e. window of 1Mb)\n\n"
	printf -- "jobs will be submitted to the cluster and flags\n"
	printf -- "--covid and --jobindex added automatically.\n\n"
	printf -- "If using a custom GWAS file, still set STUDY to\n"
	printf -- "an ID that will be used as folder and job name.\n\n"
	exit
fi

# auto-set
if [ -z "$PARQUET_FILE" ]
then
	PARQUET_FILE="${GWAS}.parquet"
fi

if [ -z "$IS_COVID" ]
then
	if [ -n "$(echo $GWAS | grep COVID)" ]
	then
		IS_COVID="--covid"
	fi
fi

# submit job
if [ -z "$JOBIND" ]
then
	printf -- "\n--------------------\n"
	printf -- "STUDY: $GWAS\n"
	printf -- "ATAC: $ATAC\n"
	printf -- "IS_COVID: $IS_COVID\n"
	printf -- "--------------------\n\n"

	N=`zcat $GEX_INFILE | wc | awk '{print $1}'`
	N=`expr '(' $N '-' 1 ')' '/' $NGENE '+' 1`

	printf "submitting jobarray with $N jobs...\n"

	mkdir -p $GWAS/LDSC
	mkdir -p $GWAS/FarmOut
	rm -f \
		$GWAS/FarmOut/stderr.%I \
		$GWAS/FarmOut/stdout.%I 
	
	MEM=20000

	temp_file=$(mktemp)
	trap "rm -f $temp_file" 0 2 3 15

	[ -n "$CUSTOM_GWAS" ] && TSV_PATH_PARAM="--gwas $CUSTOM_GWAS" || TSV_PATH_PARAM=""

	# submit job
	bsub -q normal -J"LDSCJob$GWAS[1-$N]" -e "$GWAS/FarmOut/stderr.%I" -o "$GWAS/FarmOut/stdout.%I" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" "bash getLDSC.sh $PARAMS $IS_COVID $TSV_PATH_PARAM --jobindex \$LSB_JOBINDEX" | tee $temp_file

	jobid=$(grep -Po "(?<=Job <)(\d+)(?=>)" $temp_file)

	# watch job
	watch -n 60 "bjobs -w | grep LDSCJob | grep $jobid | grep PEND | wc -l | xargs printf '"$(pwd)"\n\n> jobs for LDSCJob $GWAS $ATAC\n\nnumber of pending jobs: %s\n'; bjobs -w | grep LDSCJob | grep $jobid | grep RUN | wc -l | xargs printf 'number of running jobs: %s\n\nPress Ctrl-C to stop watching.\n'"

	exit
fi


#######################################################################
#                               RUN JOB                               #
#######################################################################

printf -- "\n--------------------\n"
printf -- "STUDY: $GWAS\n"
printf -- "ATAC: $ATAC\n"
printf -- "IS_COVID: $IS_COVID\n"
printf -- "JOBINDEX: $JOBIND\n"
printf -- "all params: $PARAMS\n"
printf -- "--------------------\n\n"

# set chunk boundaries
A=`expr '(' $JOBIND '-' 1 ')' '*' $NGENE '+' 1`
B=`expr $JOBIND '*' $NGENE`
L=`zcat $GEX_INFILE | wc | awk '{print $1}'`
if [ "$B" -ge "$L" ]
then
	B=$L
fi

# create temporary files
FLDSC="FLDSC$JOBIND.gz" #$(mktemp --suffix .gz)
FGWAS="FGWAS$JOBIND.gz" #$(mktemp --suffix .gz)

# loop over chunk
for I in `seq $A $B`
do
	TSS=`zcat $GEX_INFILE | awk -v NROW=$I 'NR==NROW {print $2}'`
	REG=`zcat $GEX_INFILE | awk -v NROW=$I -v WIN=$WSIZE '
		NR==NROW {
			if($2<WIN){
				A=WIN+1
			}else{
				A=$2
			}; 
			print $1":"A-WIN"-"$2+WIN
		}'`
	CHR=${REG%*:*}

	echo "$I, $TSS, $REG, $CHR" 1>&2

	# FLDSC
	getRsq "${VCF_FILES_DIR}chr$CHR.maf0.001.vcf.gz" \
		$REG | awk '
			BEGIN{OFS="\t"}
			$6>0.001{
				print $1"_"$2,$3,$4,$5,$6
			}' | gzip > $FLDSC

	echo "  FLDSC rows: $(wc -l $FLDSC)" 1>&2
	

	# GWAS
	if [ -n "$CUSTOM_GWAS" ]; then
		# fetch from custom file (must be tabix indexed and contain
		# columns 'hm_chrom' (1), 'hm_pos' (2+3), 'hm_other_allele' (4), 
		# 'hm_effect_allele' (5), 'hm_beta' (6) and 'standard_error' (7)
		tabix "$CUSTOM_GWAS" \
			$REG | awk -v ID=$I -v TSS=$TSS -v WIN=$WSIZE '
				BEGIN{FS="\t";OFS="\t"}
				{
					r=0.1/(0.1+$7*$7);
					bf=log(1-r)/2+($6/$7)*($6/$7)*r/2
					if($2-TSS>0){
						TSSP=($2-TSS)/WIN
					}else{
						TSSP=-($2-TSS)/WIN
					}
					print ID,$1,$2,$4,$5,bf,TSSP*(-22.91),$1"_"$2
				}' | gzip > $FGWAS
	else
		# fetch from parquet file;  old dir: /lustre/scratch117/cellgen/cellgeni/otar-gwas-ss/gwas/$GWAS.parquet; /warehouse/cellgeni/otar-gwas-ss/gwas/$GWAS.parquet
		read_parquet "$PARQUET_FILE" \
			$REG | awk -v ID=$I -v TSS=$TSS -v WIN=$WSIZE '
				BEGIN{FS="\t";OFS="\t"}
				NR>1{
					r=0.1/(0.1+$11*$11);
					bf=log(1-r)/2+($10/$11)*($10/$11)*r/2
					if($7-TSS>0){
						TSSP=($7-TSS)/WIN
					}else{
						TSSP=-($7-TSS)/WIN
					}
					print ID,$6,$7,$8,$9,bf,TSSP*(-22.91),$6"_"$7
				}' | gzip > $FGWAS
	fi

	echo "  FGWAS rows: $(wc -l $FGWAS)" 1>&2

	# join files
	join -1 8 -2 1 <(zcat $FGWAS | sort -k8,8b) <(zcat $FLDSC | sort -k1,1b) -t $'\t'
done | {
	# ATAC (optional)
	if [ -n "$ATAC" ]
	then
		n_col=$(zcat $ATAC | head -n1 | awk -F$'\t' '{printf NF; exit}')
		cellt_init=$(printf '0\t%.0s' $(seq 1 $n_col))
		cellt_init=${cellt_init%$'\t'}

		while read i || [ -n "$i" ]
		do
			LOC=${i%%$'\t'*}
			CHR=${LOC%%_*}
			POS=${LOC##*_}
			REG="$CHR:$POS-$POS"

			# query tabix indexed ATAC file
			# and summarise multiple matches
			tabix \
				$ATAC \
				$REG | awk -v LINE="$i" -v CELLT="$cellt_init" '
					BEGIN{FS="\t";OFS="\t";split(CELLT, cell_type_vec, "\t")}
					{
						# summarise in case of overlapping peak regions
						split($0, current_celltypes, "\t");
						for (i=4; i<=length(cell_type_vec); i++){
							cell_type_vec[i]=cell_type_vec[i] || current_celltypes[i];
						}
					}
					END{
						printf(LINE);
						for (i=4; i<=length(cell_type_vec); i++){
							printf(OFS cell_type_vec[i]);
						}
						printf(ORS);
					}'
		done
	else
		# if ATAC file not specified, pass input on
		cat
	fi
} | gzip > "res$JOBIND.gz"

# rm $FLDSC $FGWAS


