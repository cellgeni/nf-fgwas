# nf-fgwas
Nextflow implementation of [fGWAS](https://doi.org/10.1016/j.ajhg.2014.03.004) for single cell data. The pipeline uses the hierarchical model from [here](https://github.com/natsuhiko/PHM) and the application to single cell data has been described in detail [here](https://doi.org/10.1038/s41586-021-03852-1).

In brief, given a single cell dataset and GWAS summary statistics, fGWAS allows to find cell types that are enriched for the genetic associations.
It requires full GWAS summary statistics, including logOR and standard error.

## run

```bash
# load modules (or provide nextflow and singularity differently)
module load cellgen/nextflow  # internal module with singularity-ce version 4.1.0
module load cellgen/singularity  # internal module with nextflow version 24.10.2

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

> **Note:** if you want to use summary statistics stored on iRODS (**Sanger internal only**), run this before starting the pipeline and enter your password:
>```bash
>module load cellgen/irods
>iinit
>```

### installation

<details>
<summary><b>show details</b></summary>
<br />

`nextflow` needs to be installed to run the pipeline (obviously).
We recommend `singularity` for providing the software dependencies on an HPC environment.
A `singularity` image can be built via `docker` and a `Dockerfile` is contained in the repository.

Additionally, an already built `docker` image can be downloaded from `quay.io`:

```bash
docker pull quay.io/cellgeni/nf-fgwas
```

And then converted to `singularity`:

```bash
singularity build "nf-fgwas.sif" "docker-daemon://nf-fgwas:latest"
```

To tell nextflow to use the image, replace its path in the `nextflow.config` file.

</details>

## input files

### GWAS studies

A list of GWAS studies is supplied via a `studies.csv` file. It should contain IDs and paths to multiple summary statistics files, one study per row.

<details>
<summary><b>show details</b></summary>
<br />

Each study needs an ID (to name results folders) and summary statistics can be supplied in one of three ways:
1. path to a custom `.bed.gz` file
2. path to a `.parquet` file
3. only the study ID, in this case iRODS will be used to fetch a corresponding file called `<studyID>.parquet` from the Sanger farm (**Sanger internal only**)

The `studies.csv` file could look like this (Empty rows and lines starting with "#" will be ignored):

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

> **Note: a script is available in `bin/check_summ_stat.py` to help prepare custom files.**
See `bin/check_summ_stat.py --help` for more details.
Often summary statistics provided for different studies do not contain all required values or the columns are labelled incorrectly.
The script contains some checks to test the required values are present, as well as potentially converts them is other values are present that allow calculating the ones that are required (beta and standard error).

</details>

### RNA data

Average RNA counts per cell type need to be supplied as a TSV file. 

<details>
<summary><b>show details</b></summary>
<br />

> **Note: a script is available in `bin/prepare_input.py` to help prepare this file from an AnnData object.**
See `bin/prepare_input.py tss --help` for more details.

The TSV file should be:
- a file called `tss_cell_type_exp.txt.gz` (tab separated and gzipped) with one row per gene and:
  - 1st column: chromosome number
  - 2nd column: chromosome position
  - next columns, one for each cell type: mean expression of gene in cell type
  - last column: overall mean expression for the gene
  - sorted ascending by first two columns
  - no header

*(Currently, in addition to the gzipped file, the unzipped version including a header needs to be supplied (`--cell_types`).
This is only to provide the cell type names and will be simplified in the future.)*

</details>

### ATAC data (optional)

Optionally for including ATAC data another BED file needs to be supplied.

<details>
<summary><b>show details</b></summary>
<br />

> **Note: a script is available in `bin/prepare_input.py` to help prepare this file from an AnnData object.**
See `bin/prepare_input.py atac --help` for more details.

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

</details>

### VCF files

A path to VCF files is set in `nextflow.config`.
These are files used for calculating linkage disequilibrium scores and can be downloaded from the 1000 Genomes project website or FTP server.
(The default points to a path on Sanger HPC, replace if running from somewhere else, or other files are needed)

<details>
<summary><b>show details</b></summary>
<br />

Up-to-date (as of 2025-01-20) VCF files can be found [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/).
Fitting ancestry information: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped.

Usually, the files will have to be filtered for samples of a fitting ancestry, biallelic SNVs, and minor allele frequency >0.001.
In addition, the genome assembly (e.g. GRCh38) should also match with the summary statistics and TSS file.

Filtering the files can be done, e.g. using bcftools:
```bash
WORK_DIR="1000G_filtered_data"
VCF_URL_BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
PED_FILE_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped"
POP="EUR"  # European ancestry

mkdir -p ${WORK_DIR} && cd ${WORK_DIR}

echo "Downloading PED file for population information..."
wget -q ${PED_FILE_URL} -O 20130606_g1k.ped

echo "Extracting sample IDs for selected ancestry (${POP})..."
grep "${POP}" 20130606_g1k.ped | awk '{print $2}' > sel_samples.txt

echo "Downloading and processing VCF files..."
for CHROM in {1..22}; do
    echo "    Processing chromosome ${CHROM}..."
    VCF_FILE="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHROM}.filtered.shapeit2-duohmm-phased.vcf.gz"
    VCF_URL="${VCF_URL_BASE}/${VCF_FILE}"
    wget -q ${VCF_URL} -O ${VCF_FILE}
    
    # Filter for European samples, biallelic SNVs, and MAF > 0.001
    bcftools view -S sel_samples.txt -v snps -m2 -M2 ${VCF_FILE} | \
        bcftools view --min-af 0.001 --output-type z --output-file chr${CHROM}_filtered.vcf.gz

    # Index the filtered VCF file
    bcftools index chr${CHROM}_filtered.vcf.gz
done
```

</details>

## Citations

**this nextflow pipeline, inclusion of scATAC**

> To, Fei, Pett et al., 2024, Nature, https://doi.org/10.1038/s41586-024-08189-z

background:

**adaption for single cell data**

> Elmentaite et al., 2021, Nature, https://doi.org/10.1038/s41586-021-03852-1

**theoretical framework**

> Pickrell, 2014, AJHG, https://doi.org/10.1016/j.ajhg.2014.03.004
