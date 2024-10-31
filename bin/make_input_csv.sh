#!/bin/bash

studies="$1"
projectDir="$2"

echo "study_id,gwas_path,parquet_path" > output.csv

while IFS="" read -r var || [ -n "$var" ]
do
    if [ -z "$var" ]; then continue; fi
    if [[ ${var:0:1} != "#" ]]
    then
	study=$(echo "$var" | cut -d',' -f1)
	path=$(echo "$var" | cut -d',' -f2)

	case "$path" in
	    "$study")
		# assuming the study ID can be fetched from iRODS
		echo "$study,$projectDir/assets/NO_GWAS_FILE,$projectDir/assets/NO_PRQT_FILE" >> output.csv
		;;
	    *.parquet)
		# path to a parquet file given
		echo "$study,$projectDir/assets/NO_GWAS_FILE,$path" >> output.csv
		;;
	    *.tsv.gz)
		# path to a custom tabix indexed .tsv.gz file given
		echo "$study,$path,$projectDir/assets/NO_PRQT_FILE" >> output.csv
		;;
	    *)
		echo "Error: Unrecognized file type for study $study with path $path" >&2
		exit 1
		;;
	esac
    fi  
done < <(grep -v "^$" "${studies}")

