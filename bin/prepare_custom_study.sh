#!/bin/bash

#######################################################################
#               prepare harmonized tsv file with tabix                #
#######################################################################

FILE="$1"
OUTFILE="$(dirname $FILE)/$(basename $FILE .tsv).bed.gz"

#     |     coordinates     |                 additional fields                     |
COLS="hm_chrom hm_pos hm_pos hm_other_allele hm_effect_allele hm_beta standard_error"

# create BED file
echo "create BED file..."
cat "$FILE" | awk -v OUTCOL="$COLS" '
	BEGIN{
		FS="\t"
		OFS="\t"
		split(OUTCOL, col_names, " ")
	}
	NR==1{
		for (i=1;i<=NF;i++){
			col_pos[$i] = i
		}
	}
	NR>1{
		# exclude rows with NA values
		for (i=1; i<=length(col_names); i++){
			if ($col_pos[col_names[i]] == "NA"){
				next
			}
		}
		# keep chromosomes 1-22 (exclude X, Y and others)
		if ($col_pos[col_names[1]] ~ "^[1-9]$|^1[0-9]$|^2[0-2]$"){
			for (i=1; i<=length(col_names); i++){
				printf $col_pos[col_names[i]] (i==length(col_names)?ORS:OFS)
			}
		}
	}' | sort -t$'\t' -k1n,1 -k2n,2 | bgzip > "$OUTFILE"

# run tabix index
echo "run tabix index..."
/nfs/users/nfs_n/nk5/Applications/htslib/tabix -p bed "$OUTFILE"




