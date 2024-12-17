#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.output = true


process fetch_irods {
    input:
        val(study_id)

    output:
        tuple(val(study_id), path("${study_id}.parquet"))

    shell:
        '''
        #!/bin/bash
        archive_path=$(imeta qu -C -z archive otar_study = "!{study_id}" | sed 's/collection: //g')
        iget -fKvr "${archive_path}" .
        '''
}

process split_studies {
    input:
        path(studies)

    output:
        path('output.csv'), emit: study_data

    script:
        """
        ${projectDir}/bin/make_input_csv.sh "${studies}" "$projectDir"
        """
}

process run_LDSC {
    input:
        path("tss_cell_type_exp.txt.gz")
        path(atac_file)
        path(atac_file_tbi)
        each job_index
        val(gene_chunk_size)
        tuple(val(study_id),path(gwas_path),path(gwas_path_tbi),path(parquet_path))
  
    output:
        tuple val(study_id), path("res${job_index}.gz"), emit: res

    script:
        def use_atac_file = atac_file.name != "NO_ATAC_FILE" ? "${atac_file}" : ""
        def use_gwas_path = gwas_path.name != "NO_GWAS_FILE" ? "--gwas ${gwas_path}" : ""
        def use_parquet_path = parquet_path.name != "NO_PRQT_FILE" ? "--parquetfile ${parquet_path}" : ""
        """
        read_parquet () { "${projectDir}/bin/read_parquet.py" "\$@"; } && export -f read_parquet
        ${projectDir}/bin/getLDSC.sh "${study_id}" "${use_atac_file}" "" ${use_gwas_path} ${use_parquet_path} --jobindex "${job_index}" --vcffilesdir "${params.vcf_files_1000G}" --ngene "${gene_chunk_size}"
        """
}

process collect_LDSC {
    publishDir "results/LDSC_results", mode: "copy"

    input:
        path("tss_cell_type_exp.txt.gz")
        tuple(val(study_id), path(result_files))
        val(gene_chunk_size)

    output:
        tuple val(study_id), path("$study_id/input.gz"), emit: res

    script:
        """
        ${projectDir}/bin/collect.sh ${result_files} --without-mhc --ngene ${gene_chunk_size}
        mkdir $study_id
        mv input.gz $study_id/
        """
}

process run_HM {
    publishDir "results/HM_results", mode: "copy"

    input:
        path(atac_file)
        path(atac_file_tbi)
        path("tss_cell_type_exp.txt")
        path("tss_cell_type_exp.txt.gz")
        path("broad_fine_mapping.tsv")
        tuple(val(study_id), path("input.gz"))
        each job_index

    output:
        tuple val(study_id), path("$study_id/Out${job_index}"), emit: res

    script:
        def use_atac_file = atac_file.name != "NO_ATAC_FILE" ? "--atac ${atac_file}" : ""
        """
        ${projectDir}/bin/runHM.sh "${study_id}" ${use_atac_file} --ldscinput input.gz --jobindex "${job_index}"
        mkdir $study_id
        mv Out* $study_id/
        """
}

process plot_forest {
    publishDir "results/enrichment", mode: "copy"

    input:
        tuple(val(study_id), path(result_files))
        path("tss_cell_type_exp.txt")
        path("tss_cell_type_exp.txt.gz")

    output:
        path("$study_id/forest_plot.pdf"), emit: forest_plot
        path("$study_id/enrichment.tsv"), emit: enrichment

    script:
        """
        ${projectDir}/bin/makeForest.R --vanilla ${result_files}
        mkdir $study_id
        mv forest_plot.pdf enrichment.tsv $study_id/
        """
}

workflow fgwas {
    take:
        tss_file
        cell_types
        atac_file
        broad_fine_mapping
        gene_chunk_size
        study_data

    main:
        // ----------------- DEFINE THE INPUT CHANNELS ----------------- //

        study_data = study_data.map{ it -> tuple(it[0], it[1], "${it[1]}.tbi", it[2]) }
        atac_file_tbi = file("${atac_file}.tbi")

        // get job indices so that all rows of the tss file are processed in chunks of size gene_chunk_size
        nrow = file("${cell_types}").readLines().size() 
        log.debug "nrow: ${nrow}"
        nchunk = ((nrow - 1).intdiv(gene_chunk_size.toInteger()) + 1)
        job_indices_ngene = Channel.from(1..nchunk)

        // get job indices so that all columns (tab separated) of the tss file are processed one at a time
        ncol = file("${cell_types}").withReader{it.readLine().split("\t")}.size()
        ncelltype = (ncol - 3)
        log.debug "ncol: ${ncol}"
        job_indices_ncell = Channel.from(1..ncelltype)


        // ----------------- RUN THE WORKFLOW ----------------- //

        // 0) fetch the irods archive if no other path is given
        fetch_studies = study_data.map{ it -> 
            if (file(it[1]).name == "NO_GWAS_FILE" && file(it[3]).name == "NO_PRQT_FILE") {
                it[0]
            } else {
                null
            }
        }.filter{ it -> it != null }
        parquet_files = fetch_irods(fetch_studies)
        study_data = study_data.join(parquet_files, remainder: true).map{ it ->
            if (it[4] == null) {
                tuple(it[0], it[1], it[2], it[3])
            } else {
                tuple(it[0], it[1], it[2], it[4])
            }
        }

        // 1) run LDSC for each chunk of genes
        ldsc_results = run_LDSC(tss_file, atac_file, atac_file_tbi, job_indices_ngene, gene_chunk_size, study_data)

        // 2) collect the results of all LDSC runs
        collected_results = collect_LDSC(tss_file, ldsc_results.groupTuple(size: nchunk), gene_chunk_size)

        // 3) run HM on the collected results for all cell types
        hm_results = run_HM(atac_file, atac_file_tbi, cell_types, tss_file, broad_fine_mapping, collected_results, job_indices_ncell)

        // 4) plot the forest plot
        plot_forest(hm_results.groupTuple(size: ncelltype), cell_types, tss_file)
}

workflow {
    main:
        // ----------------- CHECK FOR REQUIRED INPUTS ----------------- //
        if (params.studies == null || !file(params.studies).exists()) {
            log.error "Missing studies file: '${params.studies}'"
            exit 1
        }
        if (!file(params.vcf_files_1000G).exists()) {
            log.error "directory not found: '${params.vcf_files_1000G}'"
            exit 1
        }
        if (params.tss_file == null || !file(params.tss_file).exists()) {
            log.error "Missing or invalid tss_file '${params.tss_file}'"
            exit 1
        }
        if (params.cell_types == null || !file(params.cell_types).exists()) {
            log.error "Missing or invalid cell_types '${params.cell_types}'"
            exit 1
        }

        // ----------------- DEFINE THE INPUT CHANNELS ----------------- //

        // define the input channels (value channels)
        tss_file = file("${params.tss_file}")  // the path to the file with transcription start sites and mean expression of genes
        cell_types = file("${params.cell_types}")  // the path to the cell type file; TODO: simplify, the file is only used for the header
        atac_file = file("${params.atac_file}")  // the path to the ATAC file (optional)
        broad_fine_mapping = file("${params.broad_fine_mapping}")  // the path to the broad fine mapping file; TODO: make optional
        gene_chunk_size = "${params.gene_chunk_size}"  // the number of genes to use per chunk / parallel job

        // split the studies file into individual studies
        study_data = split_studies(params.studies).splitCsv(header: true)
        study_data = study_data.map{it -> tuple(it.study_id, it.gwas_path, it.parquet_path)}

        study_data.tap{ study_data_tap }
        study_data_tap.view {it -> log.debug "studies: ${it}"}
        log.debug "gene chunk size: ${gene_chunk_size}"

        // ----------------- RUN THE WORKFLOW ----------------- //

        fgwas(tss_file, cell_types, atac_file, broad_fine_mapping, gene_chunk_size, study_data)
}

