#!/usr/bin/env python

import os
import logging
import tempfile
from pathlib import Path
import argparse
import numpy as np
import scipy as sp
import pandas as pd
import util
from util import add_logger


@add_logger
def fix_odds_ratio_to_beta(df, log=logging.getLogger()):
    log.info("setting hm_beta to log-transformed hm_odds_ratio.")
    df["hm_beta"] = np.sign(df["hm_odds_ratio"]) * np.log(np.abs(df["hm_odds_ratio"]))
    return df


@add_logger
def fix_harmonized_beta(df, log=logging.getLogger()):
    # TODO: better solution?
    log.warning("using beta for hm_beta. might not be harmonised!")
    df["hm_beta"] = df["beta"]
    return df


@add_logger
def fix_CI_and_beta_to_SE(df, log=logging.getLogger()):
    log.info("computing standard error from hm_beta and confidence interval.")
    df["standard_error"] = np.abs(np.log(df["hm_ci_upper"]) - np.abs(df["hm_beta"])) / 1.96
    return df


@add_logger
def fix_pvalue_and_beta_to_SE(df, log=logging.getLogger()):
    log.info("computing standard error from p-value and hm_beta.")
    df["standard_error"] = np.abs(df["hm_beta"] / sp.stats.norm.ppf(df["p_value"] / 2))
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
        'p_value': False,
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
        if col in df and not np.abs(np.nanmean(df[col]) - 0) < 0.4:
            log.error(
                f"beta values are not centered around 0, "
                f"deleting the '{col}' column."
            )
            df[col] = np.nan
            if col in required_columns:
                required_columns[col] = False
            if col in extra_columns:
                extra_columns[col] = False
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
        {"from": ["hm_odds_ratio"], "to": ["hm_beta"], "func": fix_odds_ratio_to_beta},
        {"from": ["beta"], "to": ["hm_beta"], "func": fix_harmonized_beta},
        {"from": ["hm_beta", "hm_ci_upper"], "to": ["standard_error"], "func": fix_CI_and_beta_to_SE},
        {"from": ["p_value", "hm_beta"], "to": ["standard_error"], "func": fix_pvalue_and_beta_to_SE},
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
                apply_fixes.append(fix)
                for col in fix["to"]:
                    required_columns[col] = True
                update = True  # may enable further fixes

    # apply fixes
    for fix in apply_fixes:
        df = fix["func"](df)
        for col in fix["to"]:
            if not values_present(df[col]):
                log.error(f"calculation of '{col}' unsuccessful")
                required_columns[col] = False

    
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

    try:
        # could be e.g. '10.0' or 'X' or 'NaN'
        df["hm_chrom"] = [try_conv(x) for x in df["hm_chrom"]]
        df["hm_pos"] = df.hm_pos.astype("Int64")

        log.info("reordering columns.")
        df = df[["hm_chrom", "hm_pos", "hm_pos", "hm_other_allele", "hm_effect_allele", "hm_beta", "standard_error"]]
        df.columns = ["hm_chrom", "hm_start", "hm_end", "hm_other_allele", "hm_effect_allele", "hm_beta", "standard_error"]
        log.debug(df.head())

        log.info("excluding rows with NA values in required columns.")
        df = df.dropna()
        log.debug(df.head())

        if chain_file is not None:
            # TODO: untested
            log.info(f"lifting over coordinates using chain file: '{chain_file}'")
            with tempfile.TemporaryDirectory() as tmp_dir:
                file_path = Path(tmp_dir) / "summ_stat.tsv"
                df.to_csv(file_path, sep="\t", index=False, na_rep="NA")
                
                try:
                    os.system(f'CrossMap.py bed "{chain_file}" "{file_path}" "${file_path}_tmp')
                except Error as e:
                    log.error(e)
                    log.error("is CrossMap.py installed? to install run e.g. `mamba install crossmap`")
                    raise e
                os.system(f'mv "${file_path}_tmp" "${file_path}"')
                df = pd.read_csv(file_path, sep="\t")
                log.debug(df.head())

        log.info("keeping only chromosomes 1-22.")
        df["hm_chrom"] = pd.to_numeric(df["hm_chrom"], errors='coerce', downcast='integer')
        df = df[~df["hm_chrom"].isna()]
        df["hm_chrom"] = df["hm_chrom"].astype(int)
        df = df[df["hm_chrom"].between(1, 22)]
        log.debug(df.head())

        log.info("sorting by chromosome and position.")
        df = df.sort_values(["hm_chrom", "hm_start"])
        log.debug(df.head())

    except Exception() as e:
        log.info(df.head())
        log.error(e)
        raise e

    return df


def main(log=logging.getLogger()):
    parser = argparse.ArgumentParser(description='Check and fix summary statistics obtained from EBI GWAS catalog.')

    parser.add_argument('--input_path', '-i', required=True, help='Path to the TSV file with harmonised summary statistics.')
    parser.add_argument('--output_path', '-o', default="fixed__{input_file_name}", help='Output file path for the fixed summary statistics.')
    parser.add_argument('--chain_file', default=None, help='Chain file for lifting over coordinates. No lifting over if not provided.')

    args = parser.parse_args()

    # set parameters
    output_path = Path(args.output_path.format(input_file_name=Path(args.input_path).name))
    if output_path.is_dir():
        output_path = output_path / f"fixed__{Path(args.input_path).name}"
    if output_path.suffix in [".gz", ".zip", ".bz2"]:
        output_path = output_path.with_suffix("")
    output_path = output_path.with_suffix(".tsv")
    output_path_bed = output_path.with_name(output_path.stem + ".bed.gz")
    log.info(f"output file: {output_path_bed}")

    # read the input file
    log.info("read the input file")
    df = pd.read_csv(args.input_path, sep="\t")

    # fix the summary statistics
    log.info("fix the summary statistics")
    df_fixed = fix_summ_stats(df)
    df_fixed = format_summ_stats(df_fixed, chain_file=args.chain_file)

    # save the fixed summary statistics
    log.info(f"save the fixed summary statistics: {output_path}")
    log.info(df_fixed.head())
    df_fixed.to_csv(output_path, sep="\t", index=False, na_rep="NA")

    # compress and index the output file using bgzip and tabix from htslib
    log.info("compress and index the output file using bgzip and tabix from htslib")
    with tempfile.NamedTemporaryFile() as tmp_bed:
        df_fixed.to_csv(tmp_bed, sep="\t", index=False, header=False, na_rep="NA")
        os.system(f'cat {tmp_bed.name} | bgzip > "{output_path_bed}"')
        os.system(f'tabix -p bed "{output_path_bed}"')
        log.info(f"saved tabix indexed file: {output_path_bed}")

    log.info("all done.")


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    main()
