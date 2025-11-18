#!/usr/bin/env nextflow

/*
========================================================================================
    GNPS Reference Library Annotation Pipeline
========================================================================================
    Author: Bigy Ambat
    Version: 3.0
    Description: Nextflow pipeline for GNPS reference library quality assessment
                 and spectral similarity analysis
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// Print pipeline header
def printHeader() {
    log.info """\
    ========================================
     GNPS LIBRARY ANNOTATION PIPELINE v3.0
    ========================================
    Input files    : ${params.input}
    Output dir     : ${params.outdir}

    GNPS Settings:
    m/z tolerance  : ${params.mz_tol} Da
    PPM tolerance  : ${params.ppm_tol} ppm
    Min peaks      : ${params.min_peaks}
    Cosine top N   : ${params.cosine_top_n}
    Min similarity : ${params.min_similarity}

    Resources:
    Max memory     : ${params.max_memory}
    Max CPUs       : ${params.max_cpus}

    Quality Thresholds (Research-based):
    Good score     : ${params.quality_threshold_good} (~1% FDR)
    Uncertain score: ${params.quality_threshold_uncertain} (~5% FDR)
    Min peaks (good): ${params.min_explained_peaks_good}
    Min peaks (unc.): ${params.min_explained_peaks_uncertain}
    Min peaks (min) : ${params.min_explained_peaks}
    Min intensity  : ${params.min_explained_intensity}
    ========================================
    """
}

printHeader()

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (!params.input) {
    exit 1, "Error: Please provide input GNPS MGF files using --input parameter"
}

/*
========================================================================================
    GNPS PROCESSES
========================================================================================
*/

/*
========================================================================================
    PROCESS: PARSE_GNPS_REFERENCE
========================================================================================
    Extracts reference annotations from GNPS MGF files
*/

process PARSE_GNPS_REFERENCE {
    tag "${mgf_file.baseName}"
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path mgf_file

    output:
    tuple path(mgf_file), path("${mgf_file.baseName}_reference.tsv"), emit: reference

    script:
    """
    parse_gnps_reference.py \\
        --mgf ${mgf_file} \\
        --output ${mgf_file.baseName}_reference.tsv \\
        --min_peaks ${params.min_peaks}
    """
}

/*
========================================================================================
    PROCESS: ANNOTATE_PEAKS_GNPS
========================================================================================
    Annotates peaks using GNPS reference formulas
*/

process ANNOTATE_PEAKS_GNPS {
    tag "${mgf_file.baseName}"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple path(mgf_file), path(reference_tsv)

    output:
    tuple path(mgf_file), path("${mgf_file.baseName}_annotated.tsv"), emit: annotated

    script:
    """
    annotate_peaks_gnps.py \\
        --mgf ${mgf_file} \\
        --reference ${reference_tsv} \\
        --output ${mgf_file.baseName}_annotated.tsv \\
        --mz_tolerance ${params.mz_tol} \\
        --ppm_tolerance ${params.ppm_tol} \\
        --min_explained_peaks ${params.min_explained_peaks} \\
        --min_explained_peaks_good ${params.min_explained_peaks_good} \\
        --min_explained_peaks_uncertain ${params.min_explained_peaks_uncertain} \\
        --min_explained_intensity ${params.min_explained_intensity} \\
        --quality_threshold_good ${params.quality_threshold_good} \\
        --quality_threshold_uncertain ${params.quality_threshold_uncertain}
    """
}

/*
========================================================================================
    PROCESS: CALCULATE_COSINE_SIMILARITY
========================================================================================
    Calculates pairwise cosine similarity between spectra
*/

process CALCULATE_COSINE_SIMILARITY {
    tag "${mgf_file.baseName}"
    publishDir "${params.outdir}/similarity", mode: 'copy'

    input:
    tuple path(mgf_file), path(annotations_tsv)

    output:
    tuple path(mgf_file), path("${mgf_file.baseName}_with_similarity.tsv"), emit: with_similarity

    script:
    """
    calculate_cosine_similarity.py \\
        --mgf ${mgf_file} \\
        --annotations ${annotations_tsv} \\
        --output ${mgf_file.baseName}_with_similarity.tsv \\
        --mz_tolerance ${params.mz_tol} \\
        --min_peaks ${params.min_peaks} \\
        --top_n ${params.cosine_top_n} \\
        --min_similarity ${params.min_similarity} \\
        --intensity_power 0.5
    """
}

/*
========================================================================================
    PROCESS: GENERATE_QC_GNPS
========================================================================================
    Generates HTML QC report for GNPS workflow
*/

process GENERATE_QC_GNPS {
    tag "${mgf_file.baseName}"
    publishDir "${params.outdir}/qc_reports", mode: 'copy'

    input:
    tuple path(mgf_file), path(annotations_tsv)

    output:
    path "${mgf_file.baseName}_qc_report.html", emit: qc_report

    script:
    """
    generate_qc_report_gnps.py \\
        --annotations ${annotations_tsv} \\
        --output ${mgf_file.baseName}_qc_report.html \\
        --sample_name ${mgf_file.baseName}
    """
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {
    // Create channel from input MGF files
    mgf_files_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .ifEmpty { exit 1, "No GNPS MGF files found matching pattern: ${params.input}" }

    // Display number of files to process
    mgf_files_ch
        .count()
        .subscribe { count ->
            log.info "Found ${count} GNPS MGF file(s) to process"
        }

    // Parse GNPS reference annotations
    PARSE_GNPS_REFERENCE(mgf_files_ch)

    // Annotate peaks using reference formulas
    ANNOTATE_PEAKS_GNPS(PARSE_GNPS_REFERENCE.out.reference)

    // Calculate cosine similarity
    CALCULATE_COSINE_SIMILARITY(ANNOTATE_PEAKS_GNPS.out.annotated)

    // Generate QC report
    GENERATE_QC_GNPS(CALCULATE_COSINE_SIMILARITY.out.with_similarity)

    // Log completion
    GENERATE_QC_GNPS.out.qc_report
        .collect()
        .subscribe { reports ->
            log.info """\
            ========================================
             GNPS Pipeline completed successfully!
            ========================================
             Generated ${reports.size()} QC report(s)
             Reports available in: ${params.outdir}/qc_reports
             Annotations in: ${params.outdir}/annotations
             Similarity data in: ${params.outdir}/similarity
            ========================================
            """.stripIndent()
        }
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """\
    Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Work dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline execution failed!"
    log.error "Error message: ${workflow.errorMessage}"
}
