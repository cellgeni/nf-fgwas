#!/usr/bin/env python

import logging
from collections import defaultdict
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


@add_logger
def fix_harmonize_beta(df, log=logging.getLogger()):
    # TODO: better solution?
    log.warning("rename beta to hm_beta. might not be harmonised!")
    df = df.rename(columns={"beta": "hm_beta"})


@add_logger
def fix_CI_and_beta_to_SE(df, log=logging.getLogger()):
    log.info("computing standard error from beta and confidence interval.")
    df["standard_error"] = np.abs(np.log(df["hm_ci_upper"]) - np.abs(df["hm_beta"])) / 1.96


@add_logger
def fix_pvalue_and_beta_to_SE(df, log=logging.getLogger()):
    log.info("computing standard error from p-value and beta.")
    df["standard_error"] = np.abs(df["hm_beta"] / sp.stats.norm.ppf(df["p-value"] / 2))


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
        fix_func(df)
    
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


def main():
    parser = argparse.ArgumentParser(description='Check and fix summary statistics obtained from EBI GWAS catalog.')

    parser.add_argument('--input_path', '-i', required=True, help='Path to the TSV file with harmonised summary statistics.')
    parser.add_argument('--output_path', '-o', default="fixed__{input_file_name}", help='Output file path for the fixed summary statistics (TSV)')

    args = parser.parse_args()

    # set parameters
    output_path = Path(args.output_path.format(input_file_name=Path(args.input_path).name))
    if output_path.is_dir():
        output_path = output_path / f"fixed__{Path(args.input_path).name}"

    # read the input file
    df = pd.read_csv(args.input_path, sep="\t")

    # fix the summary statistics
    df_fixed = fix_summ_stats(df)

    # save the fixed summary statistics
    df_fixed.to_csv(output_path, sep="\t", index=False)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()