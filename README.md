# Mass Spectrum Annotation & QC Pipeline

![Version](https://img.shields.io/badge/version-2.0-blue)
![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.10.0-brightgreen)
![License](https://img.shields.io/badge/license-MIT-green)

A Nextflow pipeline for automated mass spectrometry annotation and quality control reporting using MSBuddy.

## Overview

This pipeline automates the annotation and quality control reporting for mass spectrometry data. It processes `.mgf` files through MSBuddy annotation and generates comprehensive HTML QC reports with key metrics and visualizations.

### Key Features

- **Automated QC**: Eliminates manual inspection of annotation results
- **Standardized Reporting**: Consistent metrics and visualizations for every sample
- **Scalable**: Process any number of files with automatic parallelization
- **Reproducible**: Containerized dependencies ensure consistent results
- **Flexible Deployment**: Supports local, HPC cluster, and cloud execution

## Pipeline Workflow

```
MGF Files → MSBuddy Annotation → QC Report Generation → HTML Reports
```

### Processes

1. **ANNOTATE**: Runs MSBuddy on each input `.mgf` file to generate annotation results (`.tsv`)
2. **GENERATE_QC**: Creates comprehensive HTML QC reports with:
   - Overall annotation rate
   - MSBuddy score distribution
   - Precursor mass error distribution (ppm)
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
| `--mass_accuracy` | `10` | Mass accuracy in ppm |
| `--adducts` | `'[M+H]+,[M+Na]+,[M+K]+,[M-H]-,[M+Cl]-'` | Comma-separated list of adducts |
| `--ionization_mode` | `'positive'` | Ionization mode (positive/negative) |
| `--charge_range` | `'1,2'` | Charge range to consider |
| `--timeout` | `300` | Timeout per spectrum (seconds) |

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

Run with custom mass accuracy and specific adducts:

```bash
nextflow run main.nf \
    --input 'data/*.mgf' \
    --mass_accuracy 5 \
    --adducts '[M+H]+,[M+Na]+' \
    --ionization_mode positive \
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

- **Annotations** (`annotations/`): MSBuddy annotation results in TSV format
- **QC Reports** (`qc_reports/`): Self-contained HTML reports with visualizations
- **Pipeline Info** (`pipeline_info/`): Execution metrics, timeline, and DAG

## QC Report Contents

Each HTML report includes:

1. **Summary Statistics**
   - Total spectra count
   - Annotated spectra count
   - Overall annotation rate
   - Unique molecular formulas
   - Average annotation score
   - Average mass error
   - Unique adducts identified

2. **Visualizations**
   - MSBuddy score distribution histogram
   - Precursor mass error distribution histogram
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
- Version: 2.0
- Date: October 12, 2025

## Support

For issues, questions, or contributions:
- Open an issue on GitHub
- Contact: [your.email@example.com]

## Changelog

### Version 2.0 (2025-10-12)
- Initial release
- MSBuddy integration
- Automated QC report generation
- Multiple execution profiles
- Docker/Singularity support

## Future Roadmap

- [ ] MultiQC module for aggregate reporting
- [ ] Support for additional annotation engines
- [ ] Support for mzML/mzXML input formats
- [ ] Integration with metabolomics databases
- [ ] Advanced filtering and scoring options
- [ ] Publish to nf-core

## Acknowledgments

- MSBuddy development team
- Nextflow community
- Pyteomics developers
