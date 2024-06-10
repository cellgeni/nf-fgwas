#!/bin/bash

###########################
#  parquet files on farm  #
###########################

# while IFS="" read -r gwas || [ -n "$gwas" ]
# do
# 	[ -d "/lustre/scratch117/cellgen/team205/jp30/parquet_files_tmp/parquet_files/${gwas}.parquet" ] || ./copy_parquet.sh "$gwas"
# 	./submit_fGWAS.sh "$gwas" "" "" "--quiet"
# done < study_list_28Jun22_farm.txt

#########################################
#  summary stats from EBI GWAS catalog  #
#########################################

while IFS="" read -r file_path || [ -n "$file_path" ]
do
	gwas=$(echo "$file_path" | grep -Po "(?<=[^A-Z])(GCST[0-9]+)(?=[^0-9])" | head -1)
	
	echo "submit jobs for $gwas ..."

#	[ -d "EBI_downloads/${gwas}.bed.gz" ] || ./copy_EBI_summ_stat.sh "$file_path" "$gwas"

	local_file_path=$(readlink -f "$file_path")

	./submit_fGWAS.sh "$gwas" "" "$local_file_path" "--quiet"
done < study_list_8Jul22_EBI_GWAS_catalog.txt
