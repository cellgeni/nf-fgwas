#!/usr/bin/env python

import os
import logging
import tempfile
from pathlib import Path
import argparse
import numpy as np
import scipy as sp
import pandas as pd
from util import add_logger


@add_logger
def fix_odds_ratio_to_beta(df, log=logging.getLogger()):
    log.info("setting hm_beta to log-transformed hm_odds_ratio.")
    df["hm_beta"] = np.sign(df["hm_odds_ratio"]) * np.log(np.abs(df["hm_odds_ratio"]))
    return df


@add_logger
def fix_harmonize_beta(df, log=logging.getLogger()):
    # TODO: better solution?
    log.warning("rename beta to hm_beta. might not be harmonised!")
    df = df.rename(columns={"beta": "hm_beta"})
    return df


@add_logger
def fix_CI_and_beta_to_SE(df, log=logging.getLogger()):
    log.info("computing standard error from beta and confidence interval.")
    df["standard_error"] = np.abs(np.log(df["hm_ci_upper"]) - np.abs(df["hm_beta"])) / 1.96
    return df


@add_logger
def fix_pvalue_and_beta_to_SE(df, log=logging.getLogger()):
    log.info("computing standard error from p-value and beta.")
    df["standard_error"] = np.abs(df["hm_beta"] / sp.stats.norm.ppf(df["p-value"] / 2))
    return df


@add_logger
def fix_summ_stats(df, log=logging.getLogger()):
    required_columns = {
        'hm_chrom': False,
        'hm_pos': False,
        'hm_beta': False,
        'standard_error': False,
    }
    extra_columns = {
        'hm_odds_ratio': False,
        'odds_ratio': False,
        'beta': False,
        'hm_ci_upper': False,
        'p-value': False,
        'hm_ci_upper': False,
    }

    def values_present(series):
        return series.isna().mean() < 0.5
    
    # check if columns are present and have values
    for col in df.columns:
        if col in required_columns and values_present(df[col]):
            required_columns[col] = True
        if col in extra_columns and values_present(df[col]):
            extra_columns[col] = True

    # check if moving columns is necessary
    for col, to_col in {"hm_beta": "hm_odds_ratio", "beta": "odds_ratio"}.items():
        if not np.abs(np.nanmean(df[col]) - 0) < 0.4:
            log.error(
                f"beta values are not centered around 0, "
                f"deleting the '{col}' column."
            )
            del df[col]
            required_columns[col] = False
            if np.abs(np.nanmean(df[col]) - 1) < 0.4:
                log.warning(
                    "beta values are instead centered around 1, "
                    "they might be non-log-transformed OR values. "
                )
                if not required_columns[to_col]:
                    log.warning(
                        f"'{to_col}' column is missing, "
                        f"inserting the values found in '{col}'."
                    )
                    df[to_col] = df[col]
                    required_columns[to_col] = True

    # find fixes
    fixes = [
        {"from": ["hm_odds_ratio"], "to": ["hm_beta"], "func": fix_harmonize_beta},
        {"from": ["beta"], "to": ["hm_beta"], "func": fix_odds_ratio_to_beta},
        {"from": ["hm_beta", "hm_ci_upper"], "to": ["standard_error"], "func": fix_CI_and_beta_to_SE},
        {"from": ["p-value", "hm_beta"], "to": ["standard_error"], "func": fix_pvalue_and_beta_to_SE},
    ]
    update, apply_fixes = True, []
    while update:
        update = False
        for fix in fixes:
            if (
                not all(required_columns[col] for col in fix["to"]) 
                and all(required_columns[col] for col in fix["from"] if col in required_columns) 
                and all(extra_columns[col] for col in fix["from"] if col in extra_columns)
            ):
                apply_fixes.append(fix["func"])
                update = True  # may enable further fixes

    # apply fixes
    for fix_func in apply_fixes:
        df = fix_func(df)
    
    # check if all required columns are present
    if not all(required_columns.values()):
        msg = (
            "missing required columns, which could not be fixed: "
            f"{', '.join(col for col, present in required_columns.items() if not present)}"
        )
        log.error(msg)
        raise ValueError(msg)
    else:
        log.info("all done, all required columns are present.")
    
    return df


@add_logger
def format_summ_stats(df, chain_file=None, log=logging.getLogger()):
    log.info("formatting chromosome and position columns.")

    def try_conv(x):
        try:
            y=int(x)
        except ValueError:
            y=x
        return str(y)

    # could be e.g. '10.0' or 'X' or 'NaN'
    df["hm_chrom"] = [try_conv(x) for x in df["hm_chrom"]]
    df["hm_pos"] = df.hm_pos.astype("Int64")

    log.info("reordering columns.")
    df = df[["hm_chrom", "hm_pos", "hm_pos", "hm_other_allele", "hm_effect_allele", "hm_beta", "standard_error"]]

    log.info("excluding rows with NA values in required columns.")
    df = df.dropna()

    # TODO: lift over coordinates to GRCh38 if necessary (outside python, e.g. using liftOver from UCSC)
    if chain_file is not None:
        log.info(f"lifting over coordinates using chain file: '{chain_file}'")
        # save df to temporary file in temporary directory created with tempfile library
        file_path = tempfile.mktemp()
        df.to_csv(file_path, sep="\t", index=False, na_rep="NA")
        
        os.system(f'CrossMap.py bed "{chain_file}" "{file_path}" "${file_path}_tmp')
        os.system(f'mv "${file_path}_tmp" "${file_path}"')

    log.info("keeping only chromosomes 1-22.")
    df = df[df["hm_chrom"].between(1, 22)]

    log.info("sorting by chromosome and position.")
    df = df.sort_values(["hm_chrom", "hm_pos"])

    return df


def main():
    parser = argparse.ArgumentParser(description='Check and fix summary statistics obtained from EBI GWAS catalog.')

    parser.add_argument('--input_path', '-i', required=True, help='Path to the TSV file with harmonised summary statistics.')
    parser.add_argument('--output_path', '-o', default="fixed__{input_file_name}", help='Output file path for the fixed summary statistics.')
    parser.add_argument('--chain_file', default=None, help='Chain file for lifting over coordinates. No lifting over if not provided.')

    args = parser.parse_args()

    # set parameters
    output_path = Path(args.output_path.format(input_file_name=Path(args.input_path).name))
    if output_path.is_dir():
        output_path = output_path / f"fixed__{Path(args.input_path).name}"

    # read the input file
    df = pd.read_csv(args.input_path, sep="\t")

    # fix the summary statistics
    df_fixed = fix_summ_stats(df)
    df_fixed = format_summ_stats(df_fixed, chain_file=args.chain_file)

    # save the fixed summary statistics
    df_fixed.to_csv(output_path, sep="\t", index=False, na_rep="NA")

    # compress and index the output file using bgzip and tabix from htslib
    os.system(f'cat {output_path} | bgzip > "{output_path}.bed.gz"')
    os.system(f'tabix -p bed "{output_path}.bed.gz"')


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()