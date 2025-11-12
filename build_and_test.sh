#!/bin/bash

################################################################################
# GNPS Reference Library Annotation Pipeline - Build and Test Script
# Author: Bigy Ambat
# Version: 3.0
################################################################################

set -e  # Exit on error

echo "========================================="
echo " GNPS Annotation Pipeline Setup"
echo "========================================="
echo ""

# Function to display usage
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -b, --build              Build Docker image"
    echo "  -t, --test               Run test pipeline"
    echo "  -c, --clean              Clean work directory and results"
    echo "  -m, --mz-tol <value>     Set m/z tolerance in Da (default: 0.01)"
    echo "  -p, --ppm-tol <value>    Set PPM tolerance (default: 20.0)"
    echo "  -s, --min-sim <value>    Set minimum similarity (default: 0.5)"
    echo "  -h, --help               Display this help message"
    echo ""
    echo "Examples:"
    echo "  $0 --build               # Build Docker image"
    echo "  $0 --test                # Run test with Docker"
    echo "  $0 --build --test        # Build image and run test"
    echo "  $0 --test --mz-tol 0.005 --ppm-tol 5   # High-resolution data"
    echo ""
}

# Parse command line arguments
BUILD=false
TEST=false
CLEAN=false
MZ_TOL=""
PPM_TOL=""
MIN_SIM=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--build)
            BUILD=true
            shift
            ;;
        -t|--test)
            TEST=true
            shift
            ;;
        -c|--clean)
            CLEAN=true
            shift
            ;;
        -m|--mz-tol)
            MZ_TOL="$2"
            shift 2
            ;;
        -p|--ppm-tol)
            PPM_TOL="$2"
            shift 2
            ;;
        -s|--min-sim)
            MIN_SIM="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Default to help if no options provided
if [[ "$BUILD" == false && "$TEST" == false && "$CLEAN" == false ]]; then
    usage
    exit 1
fi

# Clean work directory and results
if [[ "$CLEAN" == true ]]; then
    echo "Cleaning work directory and results..."
    rm -rf work/ .nextflow/ results/
    rm -f .nextflow.log*
    echo "✓ Cleanup completed"
    echo ""
fi

# Build Docker image
if [[ "$BUILD" == true ]]; then
    echo "Building Docker image..."
    docker build -t gnps-annotation:latest .
    echo "✓ Docker image built successfully"
    echo ""
fi

# Run test
if [[ "$TEST" == true ]]; then
    echo "Checking for test data..."

    if ! ls test_data/*.mgf 1> /dev/null 2>&1; then
        echo "⚠ No GNPS MGF files found in test_data/"
        echo "Please add GNPS MGF files to test_data/ directory"
        echo "Files must contain SMILES and FORMULA fields"
        echo "See test_data/README.md for more information"
        exit 1
    fi

    echo "Found test MGF files:"
    ls -1 test_data/*.mgf
    echo ""

    # Validate MGF has required GNPS fields
    echo "Validating GNPS MGF format..."
    if ! grep -q "SMILES" test_data/*.mgf; then
        echo "⚠ Warning: No SMILES field found in MGF files"
        echo "GNPS MGF files should contain SMILES and FORMULA fields"
    fi
    if ! grep -q -E "FORMULA|MOLECULARFORMULA" test_data/*.mgf; then
        echo "⚠ Warning: No FORMULA field found in MGF files"
        echo "GNPS MGF files should contain SMILES and FORMULA fields"
    fi
    echo ""

    # Build nextflow command
    NF_CMD="nextflow run main.nf --input 'test_data/*.mgf' --outdir results -profile test,docker -resume"

    # Add custom parameters if specified
    if [[ -n "$MZ_TOL" ]]; then
        echo "Using m/z tolerance: $MZ_TOL Da"
        NF_CMD="$NF_CMD --mz_tol $MZ_TOL"
    else
        echo "Using default m/z tolerance (0.01 Da)"
    fi

    if [[ -n "$PPM_TOL" ]]; then
        echo "Using PPM tolerance: $PPM_TOL"
        NF_CMD="$NF_CMD --ppm_tol $PPM_TOL"
    else
        echo "Using default PPM tolerance (20.0)"
    fi

    if [[ -n "$MIN_SIM" ]]; then
        echo "Using minimum similarity: $MIN_SIM"
        NF_CMD="$NF_CMD --min_similarity $MIN_SIM"
    else
        echo "Using default minimum similarity (0.5)"
    fi

    echo ""
    echo "Running GNPS annotation pipeline with Docker..."
    eval $NF_CMD

    echo ""
    echo "✓ Test completed successfully"
    echo ""
    echo "Check results in:"
    echo "  - Reference:    results/reference/"
    echo "  - Annotations:  results/annotations/"
    echo "  - Similarity:   results/similarity/"
    echo "  - QC Reports:   results/qc_reports/"
    echo "  - Pipeline Info: results/pipeline_info/"
    echo ""

    # Display summary if available
    if [ -f results/qc_reports/*.html ]; then
        echo "QC Reports generated:"
        ls -1 results/qc_reports/*.html
        echo ""
        echo "Open the HTML reports in your browser to view results"
    fi
fi

echo "========================================="
echo " Setup Complete!"
echo "========================================="
