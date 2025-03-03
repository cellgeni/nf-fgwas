#!/bin/bash

studies="$1"
projectDir="$2"
launchDir="$3"

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
		case "$path" in (/*) path_tbi="$path.tbi";; (*) path_tbi="$launchDir/$path.tbi";; esac
		if [ ! -e "$path_tbi" ]; then echo "Error: no tabix index file ($path_tbi) found."; exit 1; fi
		echo "$study,$path,$projectDir/assets/NO_PRQT_FILE" >> output.csv
		;;
	    *.bed.gz)
		# path to a custom tabix indexed .bed.gz file given
		case "$path" in (/*) path_tbi="$path.tbi";; (*) path_tbi="$launchDir/$path.tbi";; esac
		if [ ! -e "$path_tbi" ]; then echo "Error: no tabix index file ($path_tbi) found."; exit 1; fi
		echo "$study,$path,$projectDir/assets/NO_PRQT_FILE" >> output.csv
		;;
	    *)
		echo "Error: Unrecognized file type for study $study with path $path" >&2
		exit 1
		;;
	esac
    fi  
done < <(grep -v "^$" "${studies}")

