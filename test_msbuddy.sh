#!/bin/bash
#
# MSBuddy Installation and Data Verification Script
# Run this on the new machine to diagnose MSBuddy issues
#

set -e

echo "========================================="
echo "MSBuddy Installation Check"
echo "========================================="

# Check Python version
echo ""
echo "1. Checking Python version..."
python --version

# Check if MSBuddy is installed
echo ""
echo "2. Checking MSBuddy installation..."
if python -c "import msbuddy" 2>/dev/null; then
    echo "   ✓ MSBuddy is installed"
    python -c "import msbuddy; print(f'   Version: {msbuddy.__version__}')"
else
    echo "   ✗ MSBuddy is NOT installed"
    echo "   Install with: pip install msbuddy"
    exit 1
fi

# Check if LightGBM is installed
echo ""
echo "3. Checking LightGBM installation..."
if python -c "import lightgbm" 2>/dev/null; then
    echo "   ✓ LightGBM is installed"
    python -c "import lightgbm; print(f'   Version: {lightgbm.__version__}')"
else
    echo "   ✗ LightGBM is NOT installed"
    echo "   Install with: pip install lightgbm"
    exit 1
fi

# Check data directory
echo ""
echo "4. Checking MSBuddy data directory..."
if [ -d "$HOME/.msbuddy" ]; then
    echo "   ✓ Data directory exists: $HOME/.msbuddy"
    echo "   Contents:"
    ls -lh "$HOME/.msbuddy/" 2>/dev/null || echo "   Directory is empty or not accessible"
else
    echo "   ⚠ Data directory does not exist: $HOME/.msbuddy"
    echo "   It may be created on first run"
fi

# Check pyteomics (needed for QC report)
echo ""
echo "5. Checking pyteomics installation..."
if python -c "import pyteomics" 2>/dev/null; then
    echo "   ✓ Pyteomics is installed"
else
    echo "   ✗ Pyteomics is NOT installed"
    echo "   Install with: pip install pyteomics"
fi

# Test MSBuddy help command
echo ""
echo "6. Testing MSBuddy command-line interface..."
if msbuddy --help > /dev/null 2>&1; then
    echo "   ✓ MSBuddy CLI is accessible"
else
    echo "   ✗ MSBuddy CLI failed"
    exit 1
fi

# Try to import and check for models
echo ""
echo "7. Checking for MSBuddy models..."
python << EOF
import msbuddy
import os
import sys

# Try to find where msbuddy is installed
msbuddy_path = msbuddy.__file__
print(f"   MSBuddy installed at: {os.path.dirname(msbuddy_path)}")

# Look for model files
model_dir = os.path.join(os.path.dirname(msbuddy_path), 'models')
if os.path.exists(model_dir):
    print(f"   ✓ Models directory found: {model_dir}")
    models = [f for f in os.listdir(model_dir) if f.endswith('.txt') or f.endswith('.model')]
    if models:
        print(f"   ✓ Found {len(models)} model file(s)")
        for m in models[:5]:  # Show first 5
            print(f"      - {m}")
    else:
        print("   ⚠ Models directory is empty")
else:
    print(f"   ⚠ Models directory not found at: {model_dir}")
    print("   Checking alternative locations...")

    # Check if models are in data directory
    data_dir = os.path.expanduser("~/.msbuddy/data")
    if os.path.exists(data_dir):
        print(f"   Found data directory: {data_dir}")
        data_files = os.listdir(data_dir)
        if data_files:
            print(f"   Contains {len(data_files)} file(s)")
EOF

echo ""
echo "========================================="
echo "Diagnostic Complete"
echo "========================================="
echo ""
echo "If all checks passed, MSBuddy should work."
echo "If there are warnings or errors, address them before running the pipeline."
