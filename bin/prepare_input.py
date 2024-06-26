#!/usr/bin/env python

### usage: `python prepare_input.py {tss,atac} --h5ad_path H5AD_PATH --groupby GROUPBY [--output_path OUTPUT_PATH]`
### help: `python prepare_input.py -h`

from pathlib import Path
import argparse
import numpy as np
import scanpy as sc
import util


def check_adata(adata, scale_thr=0, lognorm_thr=20, filter_thr=10000):
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
            adata = adata.raw.to_adata()
            adata = check_adata(adata)
        else:
            raise error
        
    return adata


def make_tss_file(h5ad_path, groupby, output_path, output_red_path, host="http://www.ensembl.org"):
    # Load the dataset
    ad = sc.read_h5ad(h5ad_path)

    # Check the dataset
    ad = check_adata(ad)
    sc.pp.filter_genes(ad, min_cells=200)

    # Compute the average expression
    avg_expr_df = util.avg_counts_by_annot(
        ad, 
        groupby, 
        add_overall_avg = "avg_expr"
    )

    # add TSS columns
    avg_expr_df = avg_expr_df.rename_axis(index="gene_name").reset_index()
    out_df = util.add_tss_to_df(
        avg_expr_df, 
        on = "gene_name",
        merge_with = "gene_symbol",
        host = host,
    )
    
    # Format and filter
    out_df = out_df.drop(columns = "gene_name")
    out_df = out_df[["gene_symbol", "gene_ensembl_id", "chromosome", "tss_loc"] + out_df.columns[4:].tolist()]
    out_df = out_df.astype({"chromosome": int, "tss_loc": int})
    out_df = out_df[out_df.chromosome.between(1, 22)]
    out_df = out_df.sort_values(["chromosome", "tss_loc"])
    out_df = out_df.drop_duplicates(["chromosome", "tss_loc"])
    
    # Save the output
    out_df.columns = util.celltypes_to_ids(out_df.columns)
    out_df.to_csv(
        Path(output_path).parent / f"gene_names_{Path(output_path).name}", 
        sep = "\t", 
        index = False,
        float_format='%.15f',
    )

    # Save the reduced output
    out_df = out_df.drop(columns = ["gene_symbol", "gene_ensembl_id"])
    out_df.to_csv(
        output_path, 
        sep = "\t", 
        index = False,
        float_format='%.15f',
    )


def make_atac_file(h5ad_path, groupby, output_path):
    # Load the dataset
    ad = sc.read_h5ad(h5ad_path)

    # Check the dataset
    ad = check_adata(ad, lognorm_thr=None)

    # Compute the accessibility
    bin_acc_df = util.get_accessibility_by_annot(
        ad, 
        groupby,
    )
    bin_acc_df["location"] = bin_acc_df["location"].str.extract("chr(.+)")[0].tolist()

    # Format and filter
    bin_acc_df = bin_acc_df.assign(
        **bin_acc_df["location"].str.extract(
            "(?P<chr>[^:\-]+):(?P<start>[^:\-]+)-(?P<end>[^:\-]+)",
            expand=True
        ).to_dict(orient="list")
    )
    bin_acc_df = bin_acc_df.iloc[:,np.r_[0, -3,-2,-1,1:-3]]
    bin_acc_df = bin_acc_df.drop(columns = "location")
    bin_acc_df = bin_acc_df.astype({"chr": int, "start": int, "end": int})
    bin_acc_df = bin_acc_df[bin_acc_df.chr.between(1, 22)]
    bin_acc_df = bin_acc_df.sort_values(["chr", "start"])

    # Save the output
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
    tss_parser.add_argument('--output_path', '-o', default="tss_cell_type_exp.txt", help='Output file path (TSV)')

    atac_parser = subparsers.add_parser('atac', help='Make ATAC file')
    atac_parser.add_argument('--h5ad_path', '-i', required=True, help='Path to the h5ad file with ATAC data')
    atac_parser.add_argument('--groupby', '-g', required=True, help='Group by column name (in AnnData.obs)')
    atac_parser.add_argument('--output_path', '-o', default="atac_cell_type_acc.bed", help='Output file path (TSV / BED)')

    args = parser.parse_args()

    if args.command == 'tss':
        make_tss_file(
            h5ad_path=args.h5ad_path,
            groupby=args.groupby,
            output_path=args.output_path,
        )
    elif args.command == 'atac':
        make_atac_file(
            h5ad_path=args.h5ad_path,
            groupby=args.groupby,
            output_path=args.output_path,
        )


if __name__ == "__main__":
    main()
