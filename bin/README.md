# run cell type enrichment

## input

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
  - list of files from open targets available on farm: `/warehouse/cellgeni/otar-gwas-ss/gwas/`
  - other summary stats can e.g. be downloaded from EBI GWAS atlas in harmonised format

optionally for including ATAC data:

- a file ending .bed.gz (tabix indexed) that contains
  - first column: chromosome position
  - second column: peak start position
  - third column: peak end position
  - remaining columns, one for each cell type:
    - 1 or 0 depending on whether peak present in cell type or not

# running

1. getLDSC.sh
2. collect.sh
3. runHM.sh
4. makeForest.sh

The script `submit_fGWAS.sh` submits all jobs in order on the farm.
