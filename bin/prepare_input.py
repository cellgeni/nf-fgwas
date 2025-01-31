#!/usr/bin/env python

### usage: `python prepare_input.py {tss,atac} --h5ad_path H5AD_PATH --groupby GROUPBY [--output_path OUTPUT_PATH]`
### help: `python prepare_input.py -h`

import sys
import os
from pathlib import Path
import logging
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import util
from util import add_logger


@add_logger
def check_adata(adata, scale_thr=0, lognorm_thr=20, filter_thr=10000, log=logging.getLogger()):
    error = None

    if scale_thr is not None and adata.X.min() < scale_thr:
        error = ValueError("adata seems to be scaled. Please provide log-transformed counts.")
    if lognorm_thr is not None and adata.X.max() > lognorm_thr:
        error = ValueError("adata seems to be raw. Please provide log-transformed counts.")
    if filter_thr is not None and adata.shape[1] < filter_thr:
        error = ValueError("adata seems to be filtered. Please provide all genes.")

    if error:
        if adata.raw is not None:
            # retry with raw attribute
            log.warning(f"AnnData check failed: '{error}'. Retrying with .raw attribute.")
            adata = adata.raw.to_adata()
            adata = check_adata(adata)
        else:
            log.error(f"AnnData check failed: '{error}' and No .raw attribute found.")
            if "Please provide log-transformed counts" in str(error):
                log.warning("log-transforming the counts.")
                sc.pp.normalize_total(adata)
                sc.pp.log1p(adata)
            else:
                raise error
        
    log.info("all checks passed")
    return adata


@add_logger
def make_tss_file(h5ad_path, groupby, output_path=None, host="http://www.ensembl.org", log=logging.getLogger()):
    # collect parameters
    if output_path is None:
        output_path = "tss_cell_type_exp.txt"
    output_path = Path(output_path)
    if Path(output_path).is_dir():
        output_path = Path(output_path) / "tss_cell_type_exp.txt"
    log.info(f"preparing TSS file...\n  h5ad file: '{h5ad_path}'\n  group by: '{groupby}'\n  output_path: '{output_path}'\n  pybiomart host: '{host}'")

    # Load the dataset
    log.info(f"load h5ad")
    ad = sc.read_h5ad(h5ad_path)

    # Check the dataset
    log.info(f"test and process")
    ad = check_adata(ad)
    sc.pp.filter_genes(ad, min_cells=200)

    # Compute the average expression
    log.info(f"average counts by group '{groupby}'")
    avg_expr_df = util.avg_counts_by_annot(
        ad, 
        groupby, 
        add_overall_avg = "avg_expr"
    )

    # add TSS columns
    log.info(f"add transcription start site info")
    avg_expr_df = avg_expr_df.rename_axis(index="gene_name").reset_index()
    out_df = util.add_tss_to_df(
        avg_expr_df, 
        on = "gene_name",
        merge_with = "gene_symbol",
        host = host,
    )
    
    # Format and filter
    log.info(f"format TSS table")
    out_df = out_df.drop(columns = "gene_name")
    out_df = out_df[["gene_symbol", "gene_ensembl_id", "chromosome", "tss_loc"] + out_df.columns[4:].tolist()]
    out_df = out_df.astype({"chromosome": int, "tss_loc": int})
    out_df = out_df[out_df.chromosome.between(1, 22)]
    out_df = out_df.sort_values(["chromosome", "tss_loc"])
    out_df = out_df.drop_duplicates(["chromosome", "tss_loc"])
    
    # Save the output with gene names (for reference)
    output_path_gene_names = Path(output_path).parent / f"gene_names_{Path(output_path).name}"
    log.info(f"save tss file with gene names and symbols for reference: '{output_path_gene_names}'")
    out_df.columns = util.celltypes_to_ids(out_df.columns)
    out_df.to_csv(
        output_path_gene_names, 
        sep = "\t", 
        index = False,
        float_format='%.15f',
    )

    # Save the output
    log.info(f"save input file for fGWAS: '{output_path}'")
    out_df = out_df.drop(columns = ["gene_symbol", "gene_ensembl_id"])
    out_df.to_csv(
        output_path, 
        sep = "\t", 
        index = False,
        float_format='%.15f',
    )


@add_logger
def make_atac_file(h5ad_path, groupby, output_path=None, metadata=None, log=logging.getLogger()):
    # collect parameters
    if output_path is None:
        output_path = "atac_cell_type_acc.bed"
    output_path = Path(output_path)
    if Path(output_path).is_dir():
        output_path = Path(output_path) / "atac_cell_type_acc.bed"
    log.info(f"preparing ATAC file...\n  h5ad file: '{h5ad_path}'\n  group by: '{groupby}'\n  output_path: '{output_path}'")

    # Load the dataset
    log.info(f"load h5ad")
    ad = sc.read_h5ad(h5ad_path)

    # Check the dataset
    log.info(f"test and process")
    ad = check_adata(ad, lognorm_thr=None)

    # add metadata if set
    if metadata is not None:
        df = pd.read_csv(metadata, index_col=0)
        ad.obs[groupby] = ad.obs_names.map(df[groupby].to_dict())

    # Compute the accessibility
    log.info(f"get accessibility by group '{groupby}'")
    bin_acc_df = util.get_accessibility_by_annot(
        ad, 
        groupby,
    )
    log.debug(bin_acc_df.head())
    bin_acc_df["location"] = bin_acc_df.index.tolist()
    bin_acc_df["location"] = bin_acc_df["location"].str.extract("chr(.+)")[0].tolist()

    # Format and filter
    log.info(f"format ATAC table")
    bin_acc_df = bin_acc_df.assign(
        **bin_acc_df["location"].str.extract(
            "(?P<chr>[^:\-]+):(?P<start>[^:\-]+)-(?P<end>[^:\-]+)",
            expand=True
        ).to_dict(orient="list")
    )
    bin_acc_df = bin_acc_df.drop(columns = "location")
    log.debug(bin_acc_df.head())

    ncol = bin_acc_df.shape[1]
    bin_acc_df = bin_acc_df.iloc[:,np.r_[ncol-3, ncol-2, ncol-1, 0:(ncol-3)]]
    log.debug(bin_acc_df.head())

    bin_acc_df.loc[~bin_acc_df["chr"].str.fullmatch("[0-9]+"), "chr"] = 0
    bin_acc_df = bin_acc_df.astype({"chr": int, "start": int, "end": int})
    bin_acc_df = bin_acc_df[bin_acc_df.chr.between(1, 22)]
    bin_acc_df = bin_acc_df.sort_values(["chr", "start"])
    log.debug(bin_acc_df.head())

    # Save the output
    log.info(f"save input file for fGWAS: '{output_path}'")
    bin_acc_df.columns = util.celltypes_to_ids(bin_acc_df.columns)
    bin_acc_df.to_csv(
        output_path, 
        sep = "\t", 
        index = False,
    )


def main():
    parser = argparse.ArgumentParser(description='Prepare input files for fGWAS pipeline')
    subparsers = parser.add_subparsers(dest='command', help='Options', title='subcommands')

    tss_parser = subparsers.add_parser('tss', help='Make TSS file')
    tss_parser.add_argument('--h5ad_path', '-i', required=True, help='Path to the h5ad file with RNA data (lognorm)')
    tss_parser.add_argument('--groupby', '-g', required=True, help='Group by column name (in AnnData.obs)')
    tss_parser.add_argument('--host', '-u', default='http://www.ensembl.org', help='ensembl URL to retrieve TSS information from. E.g. set "http://grch37.ensembl.org/" for "GRCh37.p13" or "http://jul2023.archive.ensembl.org" for "GRCh38.p14".')
    tss_parser.add_argument('--output_path', '-o', default="tss_cell_type_exp.txt", help='Output file path (TSV) or folder')
    tss_parser.set_defaults(command='tss')

    atac_parser = subparsers.add_parser('atac', help='Make ATAC file')
    atac_parser.add_argument('--h5ad_path', '-i', required=True, help='Path to the h5ad file with ATAC data')
    atac_parser.add_argument('--groupby', '-g', required=True, help='Group by column name (in AnnData.obs)')
    atac_parser.add_argument('--output_path', '-o', default="atac_cell_type_acc.bed", help='Output file path (TSV / BED) or folder')
    atac_parser.add_argument('--metadata', '-m', default=None, help='CSV file with cell metadata containing `--groupby COL` (optional)')
    atac_parser.set_defaults(command='atac')

    gzip_parser = subparsers.add_parser('gzip', help='Remove header and gzip file')
    gzip_parser.add_argument('--file', '-i', required=True, help='Path to the file (e.g. output of the tss subcommand)')
    gzip_parser.set_defaults(command='gzip')

    tabix_parser = subparsers.add_parser('tabix', help='tabix index a bed file')
    tabix_parser.add_argument('--file', '-i', required=True, help='Path to the file (e.g. output of the atac subcommand)')
    tabix_parser.set_defaults(command='tabix')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    if args.command == 'tss':
        make_tss_file(
            h5ad_path=args.h5ad_path,
            groupby=args.groupby,
            host=args.host,
            output_path=args.output_path,
        )
    elif args.command == 'atac':
        make_atac_file(
            h5ad_path=args.h5ad_path,
            groupby=args.groupby,
            output_path=args.output_path,
            metadata=args.metadata,
        )
    elif args.command == 'gzip':
        os.system(
            f"tail -n +2 {args.file} | gzip > {args.file}.gz"
        )
    elif args.command == 'tabix':
        os.system(
            f"tail -n +2 {args.file} | bgzip > {args.file}.gz && tabix -p bed {args.file}.gz"
        )


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    main()
