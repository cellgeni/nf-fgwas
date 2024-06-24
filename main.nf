#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.output = true


process fetch_irods {
    input:
        val(study_id)

    output:
        path("${study_id}.parquet"), emit: parquet_file

    shell:
        '''
        #!/bin/bash
        archive_path=$(imeta qu -C -z archive otar_study = "!{study_id}" | sed 's/collection: //g')
        iget -fKvr "${archive_path}" .
        '''
}

process run_LDSC {
    input:
        path("tss_cell_type_exp.txt.gz")
        val(study_id)
        path(atac_file)
        path(atac_file_tbi)
        path(gwas_path)
        path(gwas_path_tbi)
        path(parquet_path)
        val(job_index)
        val(gene_chunk_size)
  
    output:
        path("res${job_index}.gz"), emit: res

    publish:
        res >> null
  
    script:
        def use_atac_file = atac_file.name != "NO_ATAC_FILE" ? "${atac_file}" : ""
        def use_gwas_path = gwas_path.name != "NO_GWAS_FILE" ? "--gwas ${atac_file}" : ""
        def use_parquet_path = parquet_path.name != "NO_PRQT_FILE" ? "--parquetfile ${parquet_path}" : ""
        """
        ${projectDir}/bin/getLDSC.sh "${study_id}" "${use_atac_file}" "" ${use_gwas_path} ${use_parquet_path} --jobindex "${job_index}" --vcffilesdir "${params.vcf_files_1000G}" --ngene "${gene_chunk_size}"
        """
}

process collect_LDSC {
    input:
        path("tss_cell_type_exp.txt.gz")
        path(result_files)
        val(gene_chunk_size)

    output:
        path("input.gz"), emit: res

    publish:
        res >> "LDSC_results/"

    script:
        """
        ${projectDir}/bin/collect.sh ${result_files} --without-mhc --ngene ${gene_chunk_size}
        """
}

process run_HM {
    input:
        val(study_id)
        path(atac_file)
        path(atac_file_tbi)
        path("tss_cell_type_exp.txt")
        path("tss_cell_type_exp.txt.gz")
        path("broad_fine_mapping.tsv")
        path("input.gz")
        val(job_index)

    output:
        path("Out${job_index}"), emit: res

    publish:
        res >> "HM_results/"

    script:
        def use_atac_file = atac_file.name != "NO_ATAC_FILE" ? "--atac ${atac_file}" : ""
        """
        ${projectDir}/bin/runHM.sh "${study_id}" ${use_atac_file} --ldscinput input.gz --jobindex "${job_index}"
        """
}

process plot_forest {
    input:
        path(result_files)
        path("tss_cell_type_exp.txt")
        path("tss_cell_type_exp.txt.gz")

    output:
        path("forest_plot.pdf"), emit: forest_plot
        path("enrichment.tsv"), emit: enrichment

    publish:
        forest_plot >> "enrichment/"
        enrichment >> "enrichment/"

    script:
        """
        ${projectDir}/bin/makeForest.R --vanilla ${result_files}
        """
}

workflow {
    main:
        // ----------------- CHECK FOR REQUIRED INPUTS ----------------- //
        if (params.study_id == null) {
            log.info "Missing study_id '${params.study_id}'"
            exit 1
        }
        if (params.tss_file == null || !file(params.tss_file).exists()) {
            log.info "Missing or invalid tts_file '${params.tss_file}'"
            exit 1
        }
        if (params.cell_types == null || !file(params.cell_types).exists()) {
            log.info "Missing or invalid cell_types '${params.cell_types}'"
            exit 1
        }
        if (!file(params.vcf_files_1000G).exists()) {
            log.info "directory not found: '${params.vcf_files_1000G}'"
            exit 1
        }


        // ----------------- DEFINE THE INPUT CHANNELS ----------------- //

        // define the input channels (value channels)
        study_id = "${params.study_id}"  // the study ID (either name of a parquet file or a custom ID if gwas_path is provided)
        gwas_path = file("${params.gwas_path}")  // the path to the GWAS file (optional, may be skipped if study_id is the name of a parquet file)
        gwas_path_tbi = file("${params.gwas_path}.tbi")
        parquet_path = file("${params.parquet_path}")
        tss_file = file("${params.tss_file}")  // the path to the file with transcription start sites and mean expression of genes
        cell_types = file("${params.cell_types}")  // the path to the cell type file; TODO: simplify, the file is only used for the header
        atac_file = file("${params.atac_file}")  // the path to the ATAC file (optional)
        atac_file_tbi = file("${params.atac_file}.tbi")
        broad_fine_mapping = file("${params.broad_fine_mapping}")  // the path to the broad fine mapping file; TODO: make optional
        gene_chunk_size = "${params.gene_chunk_size}"  // the number of genes to use per chunk / parallel job

        // derive the job indices for parallel execution (queue channels)
        nrow = file("${params.cell_types}").readLines().size() 
        ncol = file("${params.cell_types}").withReader{it.readLine().split("\t")}.size()

        job_indices_ngene = Channel.from(1..((nrow - 1).intdiv(params.gene_chunk_size) + 1))  // get job indices so that all rows of the tss file are processed in chunks of size gene_chunk_size
        job_indices_ncell = Channel.from(1..ncol)  // get job indices so that all columns (tab separated) of the tss file are processed one at a time


        print tss_file
        print nrow
        print ncol


        // ----------------- RUN THE WORKFLOW ----------------- //

        // 0) fetch the irods archive
        if (gwas_path.name == "NO_GWAS_FILE" && parquet_path.name == "NO_PRQT_FILE") {
            parquet_path = fetch_irods(study_id)
        }

        // 1) run LDSC for each chunk of genes
        ldsc_results = run_LDSC(tss_file, study_id, atac_file, atac_file_tbi, gwas_path, gwas_path_tbi, parquet_path, job_indices_ngene, gene_chunk_size)

        // 2) collect the results of all LDSC runs
        collected_results = collect_LDSC(tss_file, ldsc_results.collect(), gene_chunk_size)

        // 3) run HM on the collected results for all cell types
        hm_results = run_HM(study_id, atac_file, atac_file_tbi, cell_types, tss_file, broad_fine_mapping, collected_results, job_indices_ncell)

        // 4) plot the forest plot
        plot_forest(hm_results.collect(), cell_types, tss_file)
}

output {
    directory 'results'
    enabled true
    mode 'copy'
    overwrite true
}
