# Mass Spectrum Annotation & QC Pipeline

![Version](https://img.shields.io/badge/version-2.0-blue)
![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.10.0-brightgreen)
![License](https://img.shields.io/badge/license-MIT-green)

A Nextflow pipeline for automated mass spectrometry annotation and quality control reporting using MSBuddy.

## Overview

This pipeline automates the annotation and quality control reporting for mass spectrometry data. It processes `.mgf` files through MSBuddy annotation and generates comprehensive HTML QC reports with key metrics and visualizations.

### Key Features

- **Automated QC**: Eliminates manual inspection of annotation results
- **FDR-Based Quality Classification**: Uses MSBuddy's False Discovery Rate (FDR) for reliable annotation quality assessment
- **Smart Annotation Classification**: Automatically classifies annotations as good/bad/uncertain based on statistical confidence
- **GNPS Format Support**: Properly handles GNPS library MGF files with SPECTRUMID identifiers
- **Peak Explanation Analysis**: Calculates percentage of explained peaks and intensities (when fragment data available)
- **Standardized Reporting**: Consistent metrics and visualizations for every sample
- **Configurable Quality Thresholds**: Customize criteria for good vs bad annotations
- **Scalable**: Process any number of files with automatic parallelization
- **Reproducible**: Containerized dependencies ensure consistent results
- **Flexible Deployment**: Supports local, HPC cluster, SLURM, and cloud execution

## Pipeline Workflow

```
MGF Files → MSBuddy Annotation → Peak Explanation Analysis → QC Report Generation → HTML Reports
```

### Processes

1. **ANNOTATE**: Runs MSBuddy on each input `.mgf` file to generate annotation results (`.tsv`)
2. **ANALYZE_PEAK_EXPLANATION**: Analyzes MSBuddy results to:
   - Extract FDR (False Discovery Rate) from MSBuddy predictions
   - Classify annotations as good/bad/uncertain based on FDR thresholds
   - Calculate total peaks per spectrum for quality assessment
   - Match spectrum identifiers (SPECTRUMID) between MGF and TSV files
   - Output enhanced TSV with quality metrics and classifications
3. **GENERATE_QC**: Creates comprehensive HTML QC reports with:
   - Overall annotation rate
   - Annotation quality classification breakdown
   - MSBuddy score distribution
   - Precursor mass error distribution (ppm)
   - Explained peaks distribution
   - Explained intensity distribution
   - Adduct frequency analysis

## Quick Start

### Prerequisites

- Java 8 or later (for Nextflow)
- Nextflow (>=21.10.0)
- Docker, Singularity, or Conda (for dependency management)

### Installation

1. **Install Nextflow**:
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

2. **Clone this repository**:
```bash
git clone https://github.com/yourusername/ms-annotation-qc.git
cd ms-annotation-qc
```

3. **Choose your dependency management**:

   **Option A: Docker** (Recommended)
   ```bash
   docker build -t ms-annotation-qc:latest .
   ```

   **Option B: Singularity**
   ```bash
   singularity build ms-annotation-qc.sif docker://ms-annotation-qc:latest
   ```

   **Option C: Conda**
   ```bash
   conda env create -f environment.yml
   conda activate ms-annotation-qc
   ```

### Basic Usage

```bash
# Using Docker
nextflow run main.nf \
    --input 'data/*.mgf' \
    --outdir results \
    -profile docker

# Using Singularity
nextflow run main.nf \
    --input 'data/*.mgf' \
    --outdir results \
    -profile singularity

# Using Conda
nextflow run main.nf \
    --input 'data/*.mgf' \
    --outdir results \
    -profile conda
```

## Parameters

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input` | Path pattern to input MGF files | `'data/*.mgf'` |

### Optional Parameters

#### MSBuddy Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ms1_tol` | `10` | MS1 tolerance in ppm |
| `--ms2_tol` | `10` | MS2 tolerance in ppm |
| `--timeout_secs` | `300` | Timeout per spectrum (seconds) |

#### Annotation Quality Thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--score_threshold` | `0.7` | Minimum MSBuddy score for good annotation (0-1) |
| `--explained_peaks_pct_threshold` | `0.5` | Minimum explained peaks percentage for good annotation (0-1) |
| `--explained_intensity_pct_threshold` | `0.6` | Minimum explained intensity percentage for good annotation (0-1) |
| `--mass_error_threshold` | `5.0` | Maximum mass error (ppm) for good annotation |
| `--min_explained_peaks` | `3` | Minimum number of explained peaks for good annotation |

#### Output Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `'results'` | Output directory |

#### Resource Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_memory` | `'8.GB'` | Maximum memory per process |
| `--max_cpus` | `4` | Maximum CPUs per process |
| `--max_time` | `'24.h'` | Maximum time per process |

## Examples

### Example 1: Basic Run

Process all MGF files in the current directory:

```bash
nextflow run main.nf --input '*.mgf' -profile docker
```

### Example 2: Custom Parameters

Run with custom MS tolerances and quality thresholds:

```bash
nextflow run main.nf \
    --input 'data/*.mgf' \
    --ms1_tol 5 \
    --ms2_tol 15 \
    --timeout_secs 600 \
    --score_threshold 0.8 \
    --explained_peaks_pct_threshold 0.6 \
    --mass_error_threshold 3.0 \
    --outdir my_results \
    -profile docker
```

### Example 3: HPC Cluster (SLURM)

Run on a SLURM cluster:

```bash
nextflow run main.nf \
    --input 'data/*.mgf' \
    --max_memory 16.GB \
    --max_cpus 8 \
    -profile slurm,singularity
```

### Example 4: Resume Failed Run

Resume a previously failed or interrupted run:

```bash
nextflow run main.nf --input 'data/*.mgf' -profile docker -resume
```

## Output Structure

```
results/
├── annotations/
│   ├── sample1_msbuddy.tsv
│   ├── sample2_msbuddy.tsv
│   └── ...
├── annotations_enhanced/
│   ├── sample1_enhanced.tsv
│   ├── sample2_enhanced.tsv
│   └── ...
├── qc_reports/
│   ├── sample1_qc_report.html
│   ├── sample2_qc_report.html
│   └── ...
└── pipeline_info/
    ├── execution_timeline.html
    ├── execution_report.html
    ├── execution_trace.txt
    └── pipeline_dag.svg
```

### Output Files

- **Annotations** (`annotations/`): Raw MSBuddy annotation results in TSV format
- **Enhanced Annotations** (`annotations_enhanced/`): Annotations with added quality metrics including:
  - `annotation_quality`: Classification (good/bad/uncertain/no_annotation) based on FDR
  - `estimated_fdr`: MSBuddy False Discovery Rate (primary quality metric)
  - `formula_rank_1`: Top-ranked molecular formula prediction
  - `total_peaks`: Total number of peaks in the spectrum
  - `num_explained_peaks`: Number of peaks explained (when fragment data available)
  - `explained_peaks_pct`: Percentage of explained peaks (when fragment data available)
  - `explained_intensity_pct`: Percentage of explained intensity (when fragment data available)
  - `msbuddy_score`: MSBuddy confidence score (if available)
  - `mass_error_ppm`: Mass error in ppm (if available)
- **QC Reports** (`qc_reports/`): Self-contained HTML reports with visualizations
- **Pipeline Info** (`pipeline_info/`): Execution metrics, timeline, and DAG

## Annotation Quality Classification

The pipeline uses **FDR (False Discovery Rate)** as the primary metric for classifying annotation quality. FDR represents the estimated probability that an annotation is incorrect.

### FDR-Based Classification Thresholds

| Classification | FDR Range | Meaning | Confidence Level |
|----------------|-----------|---------|------------------|
| **Good** | FDR < 0.05 | Less than 5% chance of being wrong | High confidence |
| **Uncertain** | 0.05 ≤ FDR < 0.2 | 5-20% chance of being wrong | Moderate confidence |
| **Bad** | FDR ≥ 0.2 | 20%+ chance of being wrong | Low confidence |
| **No Annotation** | N/A | No molecular formula predicted | N/A |

### Why FDR-Based Classification?

- **Statistical Rigor**: FDR provides a probabilistic measure of annotation reliability
- **MSBuddy Native**: Uses the confidence metric that MSBuddy naturally provides
- **Interpretable**: Easy to understand what FDR < 0.01 means (99% confidence)
- **Standard in Field**: FDR is widely used in metabolomics and proteomics

**Note**: MSBuddy predicts molecular formulas but does not provide MS/MS fragment explanations. Therefore, peak explanation metrics (explained peaks %, explained intensity %) are not used for quality classification with MSBuddy, but are retained for compatibility with other annotation tools.

## QC Report Contents

Each HTML report includes:

1. **Summary Statistics**
   - Total spectra count
   - Annotated spectra count
   - Overall annotation rate
   - **Good quality annotations count and percentage**
   - **Uncertain quality annotations count**
   - **Poor quality annotations count**
   - Unique molecular formulas
   - Average annotation score
   - Average mass error
   - **Average explained peaks percentage**
   - **Median explained peaks percentage**
   - **Average explained intensity percentage**
   - **Median explained intensity percentage**
   - Unique adducts identified

2. **Visualizations**
   - **Annotation quality classification bar chart**
   - MSBuddy score distribution histogram
   - Precursor mass error distribution histogram
   - **Explained peaks distribution histogram**
   - **Explained intensity distribution histogram**
   - Adduct frequency bar chart

## Execution Profiles

### Available Profiles

- `standard`: Local execution (default)
- `docker`: Docker container execution
- `singularity`: Singularity container execution
- `conda`: Conda environment execution
- `slurm`: SLURM cluster execution
- `pbs`: PBS/Torque cluster execution
- `sge`: SGE cluster execution
- `lsf`: LSF cluster execution
- `awsbatch`: AWS Batch execution
- `gcp`: Google Cloud execution
- `test`: Test profile with minimal resources
- `debug`: Debug mode with verbose output

### Profile Combination

Profiles can be combined:

```bash
nextflow run main.nf --input '*.mgf' -profile slurm,singularity
```

## Troubleshooting

### Common Issues

**Issue**: `No MGF files found matching pattern`
- **Solution**: Check your input path pattern and ensure files exist
- Use absolute paths or verify relative paths are correct

**Issue**: `MSBuddy command not found`
- **Solution**: Ensure MSBuddy is installed or use Docker/Singularity profile
- Verify PATH includes MSBuddy installation directory

**Issue**: Pipeline fails with memory errors
- **Solution**: Increase `--max_memory` parameter
- Check cluster/system resource availability

**Issue**: Reports show 'data not available' for plots
- **Solution**: Verify MSBuddy output format matches expected columns
- Check annotation TSV files are not empty
- Ensure MGF files are valid format

### Resume a Failed Pipeline

Nextflow can resume from the last successful step:

```bash
nextflow run main.nf --input '*.mgf' -resume
```

### Clean Work Directory

Remove intermediate files after successful run:

```bash
nextflow clean -f
```

## Development

### Project Structure

```
.
├── main.nf                  # Main workflow script
├── nextflow.config          # Configuration file
├── Dockerfile               # Docker image definition
├── environment.yml          # Conda environment specification
├── bin/
│   └── generate_qc_report.py  # QC report generation script
├── test_data/               # Example test data
└── README.md               # This file
```

### Testing

Run with test profile and small dataset:

```bash
nextflow run main.nf \
    --input 'test_data/*.mgf' \
    -profile test,docker
```

## Citation

If you use this pipeline in your research, please cite:

- **MSBuddy**: [MSBuddy Publication]
- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
- **pyteomics**: Goloborodko, A. A., et al. (2013). Pyteomics—a Python framework for exploratory data analysis and rapid software prototyping in proteomics. Journal of The American Society for Mass Spectrometry, 24(2), 301-304.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Author

**Bigy Ambat**
- Version: 2.1
- Date: October 30, 2025

## Support

For issues, questions, or contributions:
- Open an issue on GitHub
- Contact: [your.email@example.com]

## Changelog

### Version 2.1.1 (2025-11-04)
- **FIXED**: Spectrum ID matching to properly use SPECTRUMID field from GNPS MGF format
- **IMPROVED**: Classification logic now uses FDR (False Discovery Rate) as primary quality metric
- **NEW**: FDR-based annotation quality classification for MSBuddy:
  - `FDR < 0.01`: High confidence (good)
  - `FDR < 0.05`: Good confidence (good)
  - `FDR < 0.2`: Moderate confidence (uncertain)
  - `FDR ≥ 0.2`: Low confidence (bad)
- **NEW**: Enhanced quality statistics showing FDR distribution
- **NOTE**: MSBuddy provides molecular formula predictions but not MS/MS fragment explanations
- **NOTE**: Peak explanation metrics are retained for compatibility with other annotation tools

### Version 2.1 (2025-10-30)
- **NEW**: Peak explanation analysis with calculation of explained peaks and intensity percentages
- **NEW**: Automatic annotation quality classification (good/bad/uncertain) based on multiple criteria
- **NEW**: Configurable quality thresholds for annotation classification
- **NEW**: Enhanced TSV output with quality metrics for each spectrum
- **NEW**: Additional QC visualizations (quality classification, explained peaks, explained intensity)
- Enhanced QC reports with comprehensive quality statistics
- Updated documentation with new features

### Version 2.0 (2025-10-12)
- Initial release
- MSBuddy integration
- Automated QC report generation
- Multiple execution profiles
- Docker/Singularity support

## Future Roadmap

- [ ] MultiQC module for aggregate reporting
- [ ] GNPS library metadata integration
- [ ] Comparison with GNPS reference annotations
- [ ] Support for additional annotation engines
- [ ] Support for mzML/mzXML input formats
- [ ] Integration with metabolomics databases
- [ ] Export to SQLite database for query integration
- [ ] Publish to nf-core

## Acknowledgments

- MSBuddy development team
- Nextflow community
- Pyteomics developers
