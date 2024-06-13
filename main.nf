// TODO: publish the results somewhere
// TODO: binary used in getLDSC.sh: `/lustre/scratch126/cellgen/team205/jp30/software/vcftools/master/src/getRsq`
//       move and make sure it works on the cluster / inside singularity container
// TODO: resources used in getLDSC.sh: `/lustre/scratch126/cellgen/team205/jp30/data/genomics/1000G/20190312/European/chr$CHR.maf0.001.vcf.gz`
//       they are large files, move to a permanent location on farm and pass the path as a parameter
// TODO: binary used in getLDSC.sh: `/lustre/scratch126/cellgen/team205/jp30/software/htslib/tabix`
//       and tabix is also used through this script `/lustre/scratch126/cellgen/team205/jp30/software/scripts/tabix.py`
//       see if can replaced by one version of tabix and move to bin folder / include in environment
// TODO: binary used in runHM.sh: `/nfs/team205/jp30/projects/code/PHM/src/hm`
//       move and make sure it works on the cluster / inside singularity container

process run_LDSC {
    input:
    path("tss_cell_type_exp.txt.gz")
    val(study_id)
    path(atac_file)
    path(gwas_path)
    val(job_index)
    val(gene_chunk_size)

    output:
    path("res${job_index}.gz"), emit: res

    script:
    def use_atac_file = if (atac_file.name != "NO_FILE") "${atac_file}" else ""
    def use_gwas_path = if (gwas_path.name != "NO_FILE") "--gwas ${atac_file}" else ""
    """
    bash ${baseDir}/bin/getLDSC.sh "${study_id}" "${use_atac_file}" "" "${use_gwas_path}" --jobindex "${job_index}" --ngene "${gene_chunk_size}"
    """
}

process collect_LDSC {
    input:
    path(result_files.collect())
    val(gene_chunk_size)

    output:
    path("input.gz")

    script:
    """
    bash ${baseDir}/bin/collect.sh ${result_files} --without-mhc --ngene ${gene_chunk_size}
    """
}

process run_HM {
    input:
    val(study_id)
    path(atac_file)
    path("tss_cell_type_exp.txt")
    path("tss_cell_type_exp.txt.gz")
    path("broad_fine_mapping.tsv")
    path("input.gz")
    val(job_index)

    output:
    path("Out${job_index}"), emit: res

    script:
    def use_atac_file = if (atac_file.name != "NO_FILE") "--atac ${atac_file}" else ""
    """
    bash ${baseDir}/bin/runHM.sh "${study_id}" ${use_atac_file} --ldscinput input.gz --jobindex "${job_index}"
    """
}

process plot_forest {
    input:
    path(result_files.collect())
    path("tss_cell_type_exp.txt")
    path("tss_cell_type_exp.txt.gz")

    output:
    path("forest_plot.pdf"), emit: forest_plot
    path("enrichment.tsv"), emit: enrichment

    script:
    """
    Rscript makeForest.R --vanilla ${result_files}
    """
}

workflow {
    main:
    // ----------------- DEFINE THE INPUT CHANNELS ----------------- //

    // define the input channels (value channels)
    study_id = Channel.value("${params.study_id}")  // the study ID (either name of a parquet file or a custom ID if gwas_path is provided)
    gwas_path = path("${params.gwas_path}")  // the path to the GWAS file (optional, may be skipped if study_id is the name of a parquet file)
    tss_file = path("${params.tss_file}")  // the path to the TSS file (transcription start sites of genes and their expression levels)
    cell_types = path("${params.cell_types}")  // the path to the cell type file; TODO: simplify, the file is only used for the header
    atac_file = path("${params.atac_file}")  // the path to the ATAC file (optional)
    broad_fine_mapping = path("${params.broad_fine_mapping}")  // the path to the broad fine mapping file; TODO: make optional
    gene_chunk_size = Channel.value("${params.gene_chunk_size}")  // the number of genes to use per chunk / parallel job

    // derive the job indices for parallel execution (queue channels)
    job_indices_ngene = Channel.from(1..((tss_file.size() - 1) / gene_chunk_size + 1))  // get job indices so that all rows of the tss file are processed in chunks of size gene_chunk_size
    job_indices_ncell = Channel.from(1..tss_file.view().head().split('\t').size())  // get job indices so that all columns (tab separated) of the tss file are processed one at a time

    // ----------------- RUN THE WORKFLOW ----------------- //

    // 1) run LDSC for each chunk of genes
    ldsc_results = job_indices_ngene.map { job_index ->
        run_LDSC(tss_file, study_id, atac_file, gwas_path, job_index, gene_chunk_size)
    }

    // 2) collect the results of all LDSC runs
    collected_results = collect_LDSC(ldsc_results, gene_chunk_size)

    // 3) run HM on the collected results for all cell types
    hm_results = job_indices_ncell.map { job_index ->
        run_HM(study_id, atac_file, cell_types, tss_file, broad_fine_mapping, collected_results, job_index)
    }

    // 4) plot the forest plot
    plot_forest(hm_results, cell_types, tss_file)

    publish:
    // ----------------- PUBLISH THE RESULTS ----------------- //
    collected_results.out >> "LDSC_results/"
    hm_results.out >> "HM_results/"
    plot_forest.out >> "enrichment/"
}

output {
    directory 'results'
    enabled true
    mode 'copy'
    overwrite true
}