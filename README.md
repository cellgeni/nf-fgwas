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
nextflow run /path/to/nf-fgwas/main.nf -resume \
    --tss_file "/path/to/tss_file.txt.gz" \
    --cell_types "/path/to/tss_file.txt" \
    --atac_file "/path/to/atac_file.bed.gz" \
    --broad_fine_mapping "/path/to/broad_fine_mapping.tsv" \
    --studies "/path/to/studies.csv"
```

Alternatively, copy and edit `nextflow.config`, then execute

```bash
nextflow -C nextflow.config run /path/to/nf-fgwas/main.nf
```

**Note:** if you want to use summary statistics stored on iRODS (Sanger), also run this before starting the pipeline and enter your password:

```bash
module load cellgen/irods
iinit
```

## input files

### GWAS studies

GWAS studies are supplied via a `studies.csv` file that contains IDs and paths to multiple files with summary statistics, one study per row.
Empty rows and lines starting with "#" will be ignored. 

Each study needs an ID (to name results folders) and summary statistics can be supplied in three ways:
1. path to a custom `.bed.gz` file
2. path to a `.parquet` file
3. only the study ID, in this case iRODS will be used to fetch a file called `<studyID>.parquet` from the Sanger farm

The `studies.csv` file could look like this:

```
# fetch these studies from the farm
STUDY1
STUDY2

# parquet files
STUDY3,/path/to/study_3.parquet

# custom bed files
STUDY4,/path/to/study_4.bed.gz
```

#### custom files

Full summary statistics can often be downloaded in CSV/TSV format, e.g. from the EBI GWAS catalog in harmonised format.
Make sure that they:
   - are full summary statistics including SNP position, beta value, standard error
   - have genome build version compatible with the transcription start sites in `tss_cell_type_exp.txt.gz` (e.g. GRCh38)

Further, the files need to be transformed into BED format with these columns:
- 1st column: chromosome number (e.g. 1)
- 2nd column: SNP position (e.g. 345435)
- 3rd column: SNP position (e.g. 345435, same as 2nd column)
- 4th column: other allele (e.g. A)
- 5th column: effect allele (e.g. G)
- 6th column: beta value (e.g. 0.5)
- 7th column: standard error (e.g. 0.5)

Finally, the BED file needs to be compressed in block format using `bgzip` and indexed using `tabix`.
The index file `<file_name>.bed.gz.tbi` is expected to be in the same folder as `<file_name>.bed.gz`.

**Note: a script is available in `bin/check_summ_stat.py` to help prepare custom files.**
See `bin/check_summ_stat.py --help` for more details.
Often summary statistics provided for different studies do not contain all required values or the columns are labelled incorrectly.
The script contains some checks to test the required values are present, as well as potentially converts them is other values are present that allow calculating the ones that are required (beta and standard error).

### RNA data

Average RNA counts per cell type need to be supplied as a TSV file:

- a file called `tss_cell_type_exp.txt.gz` (tab separated and gzipped) with one row per gene and:
  - 1st column: chromosome number
  - 2nd column: chromosome position
  - next columns, one for each cell type: mean expression of gene in cell type
  - last column: overall mean expression for the gene
  - sorted ascending by first two columns
  - no header

Currently, in addition to the gzipped file, the unzipped version including a header needs to be supplied (`--cell_types`).
This is only to provide the cell type names and will be simplified in the future.

**Note: a script is available in `bin/prepare_input.py` to help prepare this file from an AnnData object.**
See `bin/prepare_input.py tss --help` for more details.

### ATAC data (optional)

Optionally for including ATAC data:

- a file ending .bed.gz (`bgzip` compressed and `tabix` indexed) that contains
  - first column: chromosome position
  - second column: peak start position
  - third column: peak end position
  - remaining columns, one for each cell type:
    - 1 or 0 depending on whether peak present in cell type or not
- a file with a cell type mapping between RNA and ATAC (`--broad_fine_mapping`)

The file provided to `--broad_fine_mapping` might look like this:
```
atac_cell_type    rna_cell_type
immune    B_cell
immune    T_cell
epithelial    Goblet
epithelial    Club
```
(note: cell type annotations can also be identical, e.g. if RNA/ATAC comes from the same cells)

**Note: a script is available in `bin/prepare_input.py` to help prepare this file from an AnnData object.**
See `bin/prepare_input.py atac --help` for more details.
