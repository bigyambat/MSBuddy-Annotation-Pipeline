#!/bin/bash

################################################################################
# GNPS Annotation Pipeline - Conda Setup (No Docker Required)
# Author: Bigy Ambat
# Version: 3.0
################################################################################

set -e  # Exit on error

echo "========================================="
echo " GNPS Pipeline - Conda Setup"
echo "========================================="
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda not found!"
    echo "Please install Miniconda or Anaconda first"
    exit 1
fi

echo "✓ Conda found: $(conda --version)"
echo ""

# Check if environment already exists
if conda env list | grep -q "^ms-annotation-qc "; then
    echo "Environment 'ms-annotation-qc' already exists"
    read -p "Do you want to recreate it? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        conda env remove -n ms-annotation-qc -y
    else
        echo "Using existing environment"
    fi
fi

# Create conda environment if it doesn't exist
if ! conda env list | grep -q "^ms-annotation-qc "; then
    echo "Creating conda environment..."
    echo "This may take 5-10 minutes..."
    conda env create -f environment.yml
    echo ""
    echo "✓ Environment created successfully"
else
    echo "✓ Environment already exists"
fi

echo ""

# Make scripts executable
echo "Making scripts executable..."
chmod +x bin/*.py
echo "✓ Scripts are now executable"
echo ""

echo "========================================="
echo " Setup Complete!"
echo "========================================="
echo ""
echo "To use the pipeline:"
echo ""
echo "1. Activate the environment:"
echo "   conda activate ms-annotation-qc"
echo ""
echo "2. Run the pipeline:"
echo "   ./run_pipeline.sh"
echo ""
echo "Or run Nextflow directly:"
echo "   nextflow run main.nf --input 'test_data/*.mgf' --outdir results"
echo ""
