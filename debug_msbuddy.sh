#!/bin/bash

################################################################################
# Debug script to test MSBuddy directly in Docker container
################################################################################

set -e

echo "========================================="
echo " MSBuddy Debug Test"
echo "========================================="
echo ""

# Check if test data exists
if ! ls test_data/*.mgf 1> /dev/null 2>&1; then
    echo "ERROR: No MGF files found in test_data/"
    exit 1
fi

TEST_MGF=$(ls test_data/*.mgf | head -1)
echo "Testing with: $TEST_MGF"
echo ""

# Test 1: Run MSBuddy directly in Docker (single-threaded, no parallel)
echo "Test 1: Running MSBuddy in Docker (single-threaded)..."
docker run --rm \
    -v $(pwd)/test_data:/data \
    -v $(pwd)/debug_output:/output \
    ms-annotation-qc:latest \
    msbuddy \
        -mgf /data/$(basename $TEST_MGF) \
        -output /output/test_output.tsv \
        -ms1_tol 10 \
        -ms2_tol 10 \
        -timeout_secs 60 \
        -n_cpu 1

echo ""
if [ -f debug_output/test_output.tsv ]; then
    echo "✓ Test 1 PASSED - Output file created"
    echo "  Lines: $(wc -l < debug_output/test_output.tsv)"
else
    echo "✗ Test 1 FAILED - No output file created"
    exit 1
fi
echo ""

# Test 2: Check MSBuddy version and help
echo "Test 2: Checking MSBuddy version..."
docker run --rm ms-annotation-qc:latest msbuddy --version || echo "Version command not available"
echo ""

# Test 3: Test with parallel mode
echo "Test 3: Running MSBuddy with parallel mode..."
docker run --rm \
    -v $(pwd)/test_data:/data \
    -v $(pwd)/debug_output:/output \
    ms-annotation-qc:latest \
    msbuddy \
        -mgf /data/$(basename $TEST_MGF) \
        -output /output/test_output_parallel.tsv \
        -ms1_tol 10 \
        -ms2_tol 10 \
        -timeout_secs 60 \
        -parallel \
        -n_cpu 2

echo ""
if [ -f debug_output/test_output_parallel.tsv ]; then
    echo "✓ Test 3 PASSED - Parallel mode works"
    echo "  Lines: $(wc -l < debug_output/test_output_parallel.tsv)"
else
    echo "✗ Test 3 FAILED - Parallel mode failed"
fi
echo ""

echo "========================================="
echo " Debug Tests Complete"
echo "========================================="
echo "Output files in: debug_output/"
