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
     MS ANNOTATION & QC PIPELINE v2.1
    ========================================
    Input files    : ${params.input}
    Output dir     : ${params.outdir}
    MS1 tolerance  : ${params.ms1_tol} ppm
    MS2 tolerance  : ${params.ms2_tol} ppm
    Timeout        : ${params.timeout_secs} seconds
    Max memory     : ${params.max_memory}
    Max CPUs       : ${params.max_cpus}

    Quality Thresholds:
    Score          : ${params.score_threshold}
    Explained peaks: ${params.explained_peaks_pct_threshold}
    Explained int. : ${params.explained_intensity_pct_threshold}
    Mass error     : ${params.mass_error_threshold} ppm
    Min peaks      : ${params.min_explained_peaks}
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
    echo "=== Starting MSBuddy annotation ==="
    echo "Input file: ${mgf_file}"
    echo "Output file: ${mgf_file.baseName}_msbuddy.tsv"
    echo "MS1 tolerance: ${params.ms1_tol} ppm"
    echo "MS2 tolerance: ${params.ms2_tol} ppm"
    echo "Timeout: ${params.timeout_secs} seconds"
    echo "CPUs: ${task.cpus}"
    echo "Memory: ${task.memory}"
    echo ""

    # Check if input file exists and is readable
    if [ ! -f "${mgf_file}" ]; then
        echo "ERROR: Input file ${mgf_file} not found!"
        exit 1
    fi

    echo "Input file size: \$(ls -lh ${mgf_file} | awk '{print \$5}')"
    echo "Input file lines: \$(wc -l < ${mgf_file})"
    echo ""

    # Run MSBuddy annotation with verbose output
    echo "Running MSBuddy..."
    msbuddy \\
        -mgf ${mgf_file} \\
        -output ${mgf_file.baseName}_msbuddy.tsv \\
        -ms1_tol ${params.ms1_tol} \\
        -ms2_tol ${params.ms2_tol} \\
        -timeout_secs ${params.timeout_secs} \\
        -parallel \\
        -n_cpu ${task.cpus}

    echo ""
    echo "=== MSBuddy annotation completed ==="
    echo "Output file size: \$(ls -lh ${mgf_file.baseName}_msbuddy.tsv | awk '{print \$5}')"
    echo "Output file lines: \$(wc -l < ${mgf_file.baseName}_msbuddy.tsv)"
    """
}

/*
========================================================================================
    PROCESS: ANALYZE_PEAK_EXPLANATION
========================================================================================
    Analyzes peak explanation and classifies annotations as good/bad
*/

process ANALYZE_PEAK_EXPLANATION {
    tag "${mgf_file.baseName}"
    publishDir "${params.outdir}/annotations_enhanced", mode: 'copy'

    input:
    tuple path(mgf_file), path(tsv_file)

    output:
    tuple path(mgf_file), path("${mgf_file.baseName}_enhanced.tsv"), emit: enhanced

    script:
    """
    analyze_peak_explanation.py \\
        --mgf ${mgf_file} \\
        --msbuddy_tsv ${tsv_file} \\
        --output ${mgf_file.baseName}_enhanced.tsv \\
        --score_threshold ${params.score_threshold} \\
        --explained_peaks_pct_threshold ${params.explained_peaks_pct_threshold} \\
        --explained_intensity_pct_threshold ${params.explained_intensity_pct_threshold} \\
        --mass_error_threshold ${params.mass_error_threshold} \\
        --min_explained_peaks ${params.min_explained_peaks}
    """
}

/*
========================================================================================
    PROCESS: GENERATE_QC
========================================================================================
    Generates HTML QC report from MGF file and enhanced annotation results
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

    // Analyze peak explanation and classify annotations
    ANALYZE_PEAK_EXPLANATION(ANNOTATE.out.annotated)

    // Run QC report generation on enhanced annotations
    GENERATE_QC(ANALYZE_PEAK_EXPLANATION.out.enhanced)

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
             Enhanced annotations in: ${params.outdir}/annotations_enhanced
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
