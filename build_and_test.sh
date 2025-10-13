#!/bin/bash

################################################################################
# MS Annotation & QC Pipeline - Build and Test Script
# Author: Bigy Ambat
# Version: 2.0
################################################################################

set -e  # Exit on error

echo "========================================="
echo " MS Annotation & QC Pipeline Setup"
echo "========================================="
echo ""

# Function to display usage
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -b, --build       Build Docker image"
    echo "  -t, --test        Run test pipeline"
    echo "  -c, --clean       Clean work directory and results"
    echo "  -h, --help        Display this help message"
    echo ""
    echo "Examples:"
    echo "  $0 --build        # Build Docker image"
    echo "  $0 --test         # Run test with Docker"
    echo "  $0 --build --test # Build image and run test"
    echo ""
}

# Parse command line arguments
BUILD=false
TEST=false
CLEAN=false

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
    docker build -t ms-annotation-qc:latest .
    echo "✓ Docker image built successfully"
    echo ""
fi

# Run test
if [[ "$TEST" == true ]]; then
    echo "Checking for test data..."

    if ! ls test_data/*.mgf 1> /dev/null 2>&1; then
        echo "⚠ No MGF files found in test_data/"
        echo "Please add test MGF files to test_data/ directory"
        echo "See test_data/README.md for more information"
        exit 1
    fi

    echo "Found test MGF files:"
    ls -1 test_data/*.mgf
    echo ""

    echo "Running test pipeline with Docker..."
    nextflow run main.nf \
        --input 'test_data/*.mgf' \
        --outdir results \
        -profile test,docker \
        -resume

    echo ""
    echo "✓ Test completed successfully"
    echo ""
    echo "Check results in:"
    echo "  - Annotations: results/annotations/"
    echo "  - QC Reports:  results/qc_reports/"
    echo "  - Pipeline Info: results/pipeline_info/"
    echo ""
fi

echo "========================================="
echo " Setup Complete!"
echo "========================================="
