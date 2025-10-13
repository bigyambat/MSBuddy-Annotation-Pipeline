#!/usr/bin/env nextflow

/*
========================================================================================
    Mass Spectrum Annotation & QC Pipeline
========================================================================================
    Author: Bigy Ambat
    Version: 2.0
    Description: Nextflow pipeline for automated MS annotation and QC reporting
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// Print pipeline header
def printHeader() {
    log.info """\
    ========================================
     MS ANNOTATION & QC PIPELINE v2.0
    ========================================
    Input files    : ${params.input}
    Output dir     : ${params.outdir}
    Mass accuracy  : ${params.mass_accuracy} ppm
    Adducts        : ${params.adducts}
    Max memory     : ${params.max_memory}
    Max CPUs       : ${params.max_cpus}
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
    exit 1, "Error: Please provide input MGF files using --input parameter"
}

/*
========================================================================================
    PROCESS: ANNOTATE
========================================================================================
    Runs MSBuddy annotation on each input MGF file
*/

process ANNOTATE {
    tag "${mgf_file.baseName}"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    path mgf_file

    output:
    tuple path(mgf_file), path("${mgf_file.baseName}_msbuddy.tsv"), emit: annotated

    script:
    """
    # Run MSBuddy annotation
    msbuddy \\
        --input ${mgf_file} \\
        --output ${mgf_file.baseName}_msbuddy.tsv \\
        --mass_accuracy ${params.mass_accuracy} \\
        --adducts ${params.adducts} \\
        --mode ${params.ionization_mode} \\
        --charge_range ${params.charge_range} \\
        --timeout ${params.timeout}
    """
}

/*
========================================================================================
    PROCESS: GENERATE_QC
========================================================================================
    Generates HTML QC report from MGF file and annotation results
*/

process GENERATE_QC {
    tag "${mgf_file.baseName}"
    publishDir "${params.outdir}/qc_reports", mode: 'copy'

    input:
    tuple path(mgf_file), path(tsv_file)

    output:
    path "${mgf_file.baseName}_qc_report.html", emit: qc_report

    script:
    """
    generate_qc_report.py \\
        --mgf ${mgf_file} \\
        --tsv ${tsv_file} \\
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
        .ifEmpty { exit 1, "No MGF files found matching pattern: ${params.input}" }

    // Display number of files to process
    mgf_files_ch
        .count()
        .subscribe { count ->
            log.info "Found ${count} MGF file(s) to process"
        }

    // Run annotation process
    ANNOTATE(mgf_files_ch)

    // Run QC report generation
    GENERATE_QC(ANNOTATE.out.annotated)

    // Log completion
    GENERATE_QC.out.qc_report
        .collect()
        .subscribe { reports ->
            log.info """\
            ========================================
             Pipeline completed successfully!
            ========================================
             Generated ${reports.size()} QC report(s)
             Reports available in: ${params.outdir}/qc_reports
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
