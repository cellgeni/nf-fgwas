# nf-fgwas
Nextflow implementation of [fGWAS](https://doi.org/10.1016/j.ajhg.2014.03.004) for single cell data. The pipeline uses the hierarchical model from [here](https://github.com/natsuhiko/PHM) and the application to single cell data has been described in detail [here](https://doi.org/10.1038/s41586-021-03852-1).

In brief, given a single cell dataset and GWAS summary statistics, fGWAS allows to find cell types that are enriched for the genetic associations.
It requires full GWAS summary statistics, including logOR and standard error.

## run

```bash
# load modules
module load cellgen/nextflow
module load cellgen/singularity

# make a new directory for the results
mkdir make/new/directory
cd make/new/directory

# run the pipeline
nextflow run /path/to/nf-fgwas/main.nf \
    --project_tag "name_of_GWAS_trait" \
    --gwas_path "/path/to/summary_statistics.bed.gz" \
    --tss_file "/path/to/tss_file.txt.gz" \
    --cell_types "/path/to/tss_file.txt" \
    --atac_file "/path/to/atac_file.bed.gz" \
    --broad_fine_mapping "/path/to/broad_fine_mapping.tsv"
```

Alternatively, copy and edit `nextflow.config`, then execute

```bash
nextflow -C nextflow.config run /path/to/nf-fgwas/main.nf
```

## requirements

- a file called `tss_cell_type_exp.txt.gz` (tab separated and gzipped) with one row per gene and:
  - 1st column: chromosome number
  - 2nd column: chromosome position
  - next columns, one for each cell type: mean expression of gene in cell type
  - last column: overall mean expression for the gene
  - sorted ascending by first two columns
  - no header

- a file called `celltype.txt`
  - list of cell type names
  - same order as in `tss_cell_type_exp.txt.gz`

- GWAS summary statistics:
  - full summary statistics including SNP position, beta value, standard error
  - genome build version compatible with the transcription start sites in `tss_cell_type_exp.txt.gz` (e.g. GRCh38)
  - list of files from open targets available on farm: `/warehouse/cellgeni/otar-gwas-ss/gwas/`
  - other summary stats can e.g. be downloaded from EBI GWAS atlas in harmonised format

optionally for including ATAC data:

- a file ending .bed.gz (tabix indexed) that contains
  - first column: chromosome position
  - second column: peak start position
  - third column: peak end position
  - remaining columns, one for each cell type:
    - 1 or 0 depending on whether peak present in cell type or not
