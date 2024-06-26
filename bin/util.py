#!/usr/bin/env python

import re
import logging
import inspect
from functools import wraps
import pandas as pd
import numpy as np
from pybiomart import Server, Dataset


def add_logger(f):
    @wraps(f)
    def wrapped(*args, **kargs):
        func_name = f.__name__

        log = logging.getLogger(func_name)
        if not log.hasHandlers():
            c_handler = logging.StreamHandler()
            c_handler.setLevel(logging.INFO)
            c_format = logging.Formatter(
                '[%(levelname)s - %(name)s - %(asctime)s]: %(message)s'
            )
            c_handler.setFormatter(c_format)
            log.addHandler(c_handler)

        log.info(f"----- starting {func_name} -----")
        r = f(*args, **kargs, log=log)
        log.info("----- all done -----")

        return r
    # edit signature to remove params (PEP-0362)
    remove_params = ["log"]
    sig = inspect.signature(f)
    sig = sig.replace(parameters=tuple(v for k,v in sig.parameters.items() if k not in remove_params))
    wrapped.__signature__ = sig
    return wrapped


@add_logger
def celltypes_to_ids(
    cell_type_list,
    log = logging.getLogger(),
):
    id_map = {
        "\+": "pos",
        "/": "or",
        "[^A-Za-z0-9_]+": "_",
    }
    
    log.debug(id_map)
    
    cellt_ids = []
    
    for ct in cell_type_list:
        ct_id = ct
        for m, r in id_map.items():
            ct_id = re.sub(m, r, ct_id)
        cellt_ids.append(ct_id)
        
    return cellt_ids


@add_logger
def avg_counts_by_annot(
    ad, 
    annot,
    add_overall_avg = None,
    log = logging.getLogger(),
):
    """
    Average counts in anndata `ad.X` per annotation `annot`, 
    where `annot` is a column name of `ad.obs`.
    Return pandas DataFrame (genes x annot)
    """

    annot_cat = ad.obs[annot].unique().tolist()

    avg_expr_df = pd.DataFrame(index = ad.var_names, columns = annot_cat)

    for ac in annot_cat:
        log.info(f"averaging counts for {ac}")
        mask = ad.obs[annot] == ac
        avg_counts = ad.X[mask,:].mean(axis=0)
        avg_expr_df[ac] = np.array(avg_counts).flatten()
        
    if add_overall_avg:
        avg_expr_df[add_overall_avg] = avg_expr_df.mean(axis=1)
        
    return avg_expr_df


@add_logger
def get_accessibility_by_annot(
    ad, 
    annot,
    frac_thr = 0.1,
    log = logging.getLogger(),
):
    """
    Get binary accessibility based on a threshold `frac_thr` on 
    the fraction of accessible cells for anndata `ad.X` per 
    annotation `annot`, where `annot` is a column name of `ad.obs`.
    Return pandas DataFrame (peaks x annot)
    """

    annot_cat = ad.obs[annot].unique().tolist()

    bin_acc_df = pd.DataFrame(index = ad.var_names, columns = annot_cat)

    for ac in annot_cat:
        log.info(f"getting accessibility for {ac}")
        mask = ad.obs[annot] == ac
        n_cells = mask.sum()
        n_cells_acc = (ad.X[mask,:] > 0).sum(axis=0)
        
        acc_after_thr = (n_cells_acc / n_cells) > frac_thr
        
        bin_acc_df[ac] = np.array(acc_after_thr).flatten().astype(int)
        
    return bin_acc_df


def get_tss_df(host=None):
    host = host or 'http://www.ensembl.org'
    # load biomart dataset
    server = Server(host=host)
    mart = server['ENSEMBL_MART_ENSEMBL']
    dataset = mart['hsapiens_gene_ensembl']
    
    # extract data frame
    def select_tss(df):
        strand = df["Strand"].value_counts().index[0]
        chromosome = df["Chromosome/scaffold name"].value_counts().index[0]
        gene_symbol = df["Gene name"].value_counts().index[0]
        gene_ensembl_id = df["Gene stable ID"].value_counts().index[0]
        tss = df["Transcription start site (TSS)"].agg(min if strand > 0 else max)
        return pd.DataFrame({
            "chromosome": [chromosome], 
            "tss_loc": [tss],
            "gene_symbol": [gene_symbol],
            "gene_ensembl_id": [gene_ensembl_id],
        })
    
    tss_df = dataset.query(
        attributes = [
            'transcription_start_site',
            'strand',
            'external_gene_name',
            'ensembl_gene_id',
            'chromosome_name'
        ]
    ).groupby("Gene name").apply(select_tss).droplevel(1)
    
    return tss_df


def _check_opt(val, arg):
    opt_dict = {
        "merge_with": ["gene_symbol", "gene_ensembl_id"],
    }
    if val not in opt_dict[arg]: 
        opt_str = ', '.join(f"'{x}'" for x in opt_dict[arg])
        raise ValueError(f"'{arg}' options: [{opt_str}]")


def merge_tss_into_df(df, tss_df, on="gene_symbol", merge_with="gene_symbol"):
    _check_opt(merge_with, 'merge_with')
    
    mrg_df = tss_df.merge(
        df, 
        how = 'inner',
        left_on = merge_with,
        right_on = on,
    ).query(
        "chromosome.str.match('[0-9]+$')", 
        engine='python'
    ).astype({
        "chromosome": int, 
        "tss_loc": int,
    }).sort_values(
        ["chromosome", "tss_loc"]
    )
    
    return mrg_df


@add_logger
def add_tss_to_df(
        df, 
        on="gene_symbol", 
        merge_with="gene_symbol", 
        host='http://www.ensembl.org',
        log = logging.getLogger(),
    ):
    _check_opt(merge_with, 'merge_with')

    log.info(f"getting from host: {host}")
    tss_df = get_tss_df(host=host)
    
    log.info("adding TSS to DataFrame")
    mrg_df = merge_tss_into_df(
        df, 
        tss_df, 
        on = on, 
        merge_with = merge_with,
    )
    
    return mrg_df