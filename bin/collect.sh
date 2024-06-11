#!/bin/bash

#######################################################################
#                          HANDLE PARAMETERS                          #
#######################################################################

# fixed parameters
GEX_INFILE="tss_cell_type_exp.txt.gz"
OUTFILE="input.gz"

# default parameters
NGENE=100
WOMHC=""
PARAMS=""

while (( "$#" ))
do
	case "$1" in
		-w|--without-mhc)
			# skip MHC region
			WOMHC=true
			shift
			;;
		-n|--ngene)
			# set chunk size
			if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
				NGENE="$2"
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

FILE_CHUNKS=("${PARGS[@]}")  # all positional arguments are files to be merged


#######################################################################
#                             RUN SCRIPT                              #
#######################################################################

# get number of chunks
N=`zcat $GEX_INFILE | wc | awk '{print $1}'`
N=`expr '(' $N '-' 1 ')' '/' $NGENE '+' 1`


echo "starting to collect files..." >&2

I=1
for p in "${FILE_CHUNKS[@]}"
	do
		echo "$I / $N --> $p" >&2
		I=$((I + 1))
		zcat "$p" | awk '
			BEGIN{OFS="\t"}
			$12>0.001{
				printf "%s%s%s%s%s", $2, OFS, $7, OFS, $8-log($12);
				for (i=13;i<=NF;i++) printf "%s%s", OFS, $i;
				printf "%s", ORS;
			}'
	done | gzip > "$OUTFILE"


if [ -n "$WOMHC" ]
then
	# TODO: these regions are hard coded, should be passed as parameters
	MHCA=`zcat $GEX_INFILE | awk '$1==6&&$2>28510120&&$2<33480577{print NR}' | head -n 1`
	MHCB=`zcat $GEX_INFILE | awk '$1==6&&$2>28510120&&$2<33480577{print NR}' | tail -n 1`

	zcat "$OUTFILE" | awk -v MHCA=$MHCA -v MHCB=$MHCB '$1<MHCA || $1>MHCB' | gzip > "$OUTFILE"
fi


zcat "$OUTFILE" | wc | awk '{print $1}'


#| awk 'BEGIN{I=0;J=0;OFS="\t"}{if($1!=I){I=$1;J=J+1}else{};print J,$1,$2,$3}' | gzip > input.gz

