#!/bin/bash

# input gex_infile: <tmp/tss_cell_type_exp.txt.gz>
# output collected_chunks_out: <getLDSC/collected/input_hm.gz>

GEX_INFILE="tmp/tss_cell_type_exp.txt.gz"

if [ -z $1 ]
then
	echo Specify study!
	exit
fi

CHUNK_OUT="$3"
NGENE="$4"



MHCA=`zcat $GEX_INFILE | awk '$1==6&&$2>28510120&&$2<33480577{print NR}' | head -n 1`
MHCB=`zcat $GEX_INFILE | awk '$1==6&&$2>28510120&&$2<33480577{print NR}' | tail -n 1`

N=`zcat $GEX_INFILE | wc | awk '{print $1}'`
N=`expr '(' $N '-' 1 ')' '/' $NGENE '+' 1`

echo "starting to collect files..." >&2

if [ -z $2 ]
then
	for I in `seq 1 $N`
	do
		echo "$I / $N --> getLDSC/chunks/res$I.gz" >&2
		zcat getLDSC/chunks/res$I.gz | awk '
			BEGIN{OFS="\t"}
			$12>0.001{
				printf "%s%s%s%s%s", $2, OFS, $7, OFS, $8-log($12);
				for (i=13;i<=NF;i++) printf "%s%s", OFS, $i;
				printf "%s", ORS;
			}'
	done | gzip > "$CHUNK_OUT"
	zcat "$CHUNK_OUT" | wc | awk '{print $1}'
else
	for I in `seq 1 $N`
	do
		echo "$I / $N --> getLDSC/chunks/res$I.gz" >&2
		zcat getLDSC/chunks/res$I.gz | awk '
			BEGIN{OFS="\t"}
			$12>0.001{
				printf "%s%s%s%s%s", $2, OFS, $7, OFS, $8-log($12);
				for (i=13;i<=NF;i++) printf "%s%s", OFS, $i;
				printf "%s", ORS;
			}'
 
	done | awk -v MHCA=$MHCA -v MHCB=$MHCB '$1<MHCA || $1>MHCB' | gzip > "$CHUNK_OUT"
	zcat "$CHUNK_OUT" | wc | awk '{print $1}'
fi

#| awk 'BEGIN{I=0;J=0;OFS="\t"}{if($1!=I){I=$1;J=J+1}else{};print J,$1,$2,$3}' | gzip > input.gz

