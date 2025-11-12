# GNPS Reference Library Annotation Pipeline

![Version](https://img.shields.io/badge/version-3.0-blue)
![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.10.0-brightgreen)
![License](https://img.shields.io/badge/license-MIT-green)

A Nextflow pipeline for GNPS reference library quality assessment and spectral similarity analysis.

## Overview

This pipeline processes GNPS reference library MGF files to:

1. **Extract reference annotations** (SMILES, formulas) from MGF metadata
2. **Annotate fragment peaks** using reference molecular formulas
3. **Classify annotation quality** (good/uncertain/bad) based on peak explanation
4. **Calculate spectral similarity** using cosine scores
5. **Generate interactive QC reports** with comprehensive visualizations

### Key Features

- **Quality Assessment**: Separates high-quality from low-quality reference annotations
- **Spectral Clustering**: Identifies similar spectra for molecular networking
- **Cosine Similarity**: Finds spectrally similar neighbors for metadata imputation
- **Interactive Reports**: Self-contained HTML with 8 visualization types
- **Scalable**: Handles datasets from 100 to 100,000+ spectra
- **Flexible Deployment**: Local, HPC cluster, Docker, Singularity, or cloud execution

## Pipeline Workflow

```
GNPS MGF Files (with SMILES & Formulas)
    ↓
Parse Reference Annotations
    ↓
Annotate Peaks → Calculate Quality Metrics
    ↓
Calculate Cosine Similarity
    ↓
Generate Interactive HTML QC Report
```

### Processes

1. **PARSE_GNPS_REFERENCE**: Extracts SMILES, InChI, formulas from GNPS MGF metadata
2. **ANNOTATE_PEAKS_GNPS**: Matches theoretical fragments to observed peaks, calculates quality scores
3. **CALCULATE_COSINE_SIMILARITY**: Computes pairwise spectral similarities, finds top neighbors
4. **GENERATE_QC_GNPS**: Creates comprehensive HTML reports with visualizations

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
git clone https://github.com/yourusername/gnps-annotation-pipeline.git
cd gnps-annotation-pipeline
```

3. **Choose your dependency management**:

   **Option A: Docker** (Recommended)
   ```bash
   docker build -t gnps-annotation:latest .
   ```

   **Option B: Singularity**
   ```bash
   singularity build gnps-annotation.sif docker://gnps-annotation:latest
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
    --input 'gnps_data/*.mgf' \
    --outdir results \
    -profile docker

# Using Singularity
nextflow run main.nf \
    --input 'gnps_data/*.mgf' \
    --outdir results \
    -profile singularity

# Using Conda
nextflow run main.nf \
    --input 'gnps_data/*.mgf' \
    --outdir results \
    -profile conda
```

## Input Format

### GNPS MGF Requirements

Your MGF files must contain:

```
BEGIN IONS
TITLE=Spectrum_Name
SPECTRUMID=CCMSLIB00000001234        ← GNPS library ID
SMILES=CC(C)CC1=CC=C(C=C1)C(C)C      ← Required: SMILES structure
FORMULA=C15H24                        ← Required: Molecular formula
PEPMASS=205.1956                      ← Required: Precursor m/z
CHARGE=1+
123.4567 1000.0                       ← Fragment peaks (m/z intensity)
234.5678 500.0
END IONS
```

**Required Fields:**
- `SMILES`: Chemical structure
- `FORMULA` or `MOLECULARFORMULA`: Molecular formula
- `PEPMASS`: Precursor m/z
- Fragment peak list

**Optional Fields:**
- `SPECTRUMID`, `INCHI`, `INCHIKEY`, `NAME`, `COLLISION_ENERGY`, `RTINSECONDS`

## Parameters

### Peak Matching

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mz_tol` | `0.01` | m/z tolerance for peak matching (Da) |
| `--ppm_tol` | `20.0` | PPM tolerance for peak matching |
| `--min_peaks` | `3` | Minimum peaks required in spectrum |

### Cosine Similarity

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--cosine_top_n` | `10` | Number of top similar neighbors to report |
| `--min_similarity` | `0.5` | Minimum cosine similarity score (0-1) |
| `--intensity_power` | `0.5` | Power scaling for intensity (0.5 = sqrt) |

### Quality Thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--quality_threshold_good` | `0.7` | Minimum score for "good" classification |
| `--quality_threshold_uncertain` | `0.4` | Minimum score for "uncertain" |
| `--min_explained_peaks` | `3` | Minimum explained peaks required |
| `--min_explained_intensity` | `0.5` | Minimum explained intensity fraction |

### Resource Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_memory` | `16.GB` | Maximum memory per process |
| `--max_cpus` | `8` | Maximum CPUs per process |
| `--max_time` | `24.h` | Maximum runtime per process |

## Output Structure

```
results/
├── reference/                        ← Parsed reference annotations
│   └── sample_reference.tsv
├── annotations/                      ← Peak annotations with quality metrics
│   └── sample_annotated.tsv
├── similarity/                       ← Cosine similarity results
│   ├── sample_with_similarity.tsv   (Enhanced annotations)
│   └── sample_neighbors.tsv         (Top N neighbors per spectrum)
├── qc_reports/                       ← Interactive HTML reports
│   └── sample_qc_report.html
└── pipeline_info/                    ← Nextflow execution logs
    ├── execution_timeline.html
    ├── execution_report.html
    └── execution_trace.txt
```

### Key Output Files

**1. Reference Annotations** (`reference/{sample}_reference.tsv`)
- Parsed SMILES, InChI, formulas, metadata
- Validation flags: `has_structure`, `has_formula`, `annotation_complete`

**2. Annotated Peaks** (`annotations/{sample}_annotated.tsv`)
- All reference data plus:
- `num_matched_peaks`, `matched_peaks_pct`, `matched_intensity_pct`
- `quality_score` (0-1), `annotation_quality` (good/uncertain/bad)

**3. Similarity-Enhanced** (`similarity/{sample}_with_similarity.tsv`)
- All annotation data plus:
- `num_similar_neighbors`, `max_similarity_score`, `avg_similarity_score`
- `top_neighbor`, `top_neighbor_score`

**4. Top Neighbors** (`similarity/{sample}_neighbors.tsv`)
- Pairwise similarity relationships
- `query_spectrum`, `match_spectrum`, `rank`, `cosine_score`

**5. QC Report** (`qc_reports/{sample}_qc_report.html`)
- Interactive HTML with 8 plots:
  - Annotation quality distribution
  - Quality score distribution
  - Matched peaks/intensity distributions
  - Spectral similarity distribution
  - Quality vs. similarity scatter plot
  - Similar neighbors distribution
  - Adduct frequency

## Usage Examples

### High-Resolution Data (Orbitrap, FTICR)

```bash
nextflow run main.nf \
    --input 'orbitrap_data/*.mgf' \
    --outdir results \
    --mz_tol 0.005 \
    --ppm_tol 5.0 \
    -profile docker
```

### Low-Resolution Data (Ion Trap, TOF)

```bash
nextflow run main.nf \
    --input 'qtof_data/*.mgf' \
    --outdir results \
    --mz_tol 0.05 \
    --ppm_tol 50.0 \
    -profile docker
```

### Large Dataset (HPC Cluster)

```bash
nextflow run main.nf \
    --input 'large_library/*.mgf' \
    --outdir results \
    --min_similarity 0.7 \
    --max_memory 32.GB \
    -profile slurm,singularity
```

### Strict Quality Filtering

```bash
nextflow run main.nf \
    --input 'gnps_data/*.mgf' \
    --outdir results \
    --quality_threshold_good 0.8 \
    --min_explained_peaks 5 \
    --min_explained_intensity 0.7 \
    -profile docker
```

## Advanced Analysis

### Filter High-Quality Annotations

```python
import pandas as pd

# Load results
df = pd.read_csv('results/similarity/sample_with_similarity.tsv', sep='\t')

# Filter good annotations
good_annotations = df[
    (df['annotation_quality'] == 'good') &
    (df['matched_intensity_pct'] > 70) &
    (df['smiles'].notna())
]

# Save filtered set
good_annotations.to_csv('high_quality_library.tsv', sep='\t', index=False)

print(f"High quality: {len(good_annotations)} / {len(df)} ({len(good_annotations)/len(df)*100:.1f}%)")
```

### Build Spectral Network

```python
import pandas as pd
import networkx as nx

# Load similarity neighbors
neighbors = pd.read_csv('results/similarity/sample_neighbors.tsv', sep='\t')

# Build network (edges = spectral similarity)
G = nx.from_pandas_edgelist(
    neighbors,
    source='query_spectrum',
    target='match_spectrum',
    edge_attr='cosine_score'
)

# Find clusters (connected components)
clusters = list(nx.connected_components(G))
print(f"Found {len(clusters)} spectral clusters")

# Analyze largest cluster
largest_cluster = max(clusters, key=len)
print(f"Largest cluster: {len(largest_cluster)} spectra")
```

### Impute Missing Metadata

```python
import pandas as pd
import numpy as np

annotations = pd.read_csv('results/similarity/sample_with_similarity.tsv', sep='\t')
neighbors = pd.read_csv('results/similarity/sample_neighbors.tsv', sep='\t')

# Impute collision energy from similar spectra
for idx, row in annotations[annotations['collision_energy'].isna()].iterrows():
    spec_id = row['spectrum_id']

    # Get similar neighbors
    similar = neighbors[neighbors['query_spectrum'] == spec_id]
    neighbor_ids = similar['match_spectrum'].tolist()

    # Find neighbors with known CE
    neighbor_data = annotations[
        (annotations['spectrum_id'].isin(neighbor_ids)) &
        (annotations['collision_energy'].notna())
    ]

    if len(neighbor_data) > 0:
        # Weighted average by similarity score
        weights = similar[similar['match_spectrum'].isin(neighbor_data['spectrum_id'])]['cosine_score']
        imputed_ce = np.average(neighbor_data['collision_energy'], weights=weights)
        annotations.at[idx, 'collision_energy'] = imputed_ce
        annotations.at[idx, 'ce_imputed'] = True

# Save enhanced annotations
annotations.to_csv('annotations_with_imputed_metadata.tsv', sep='\t', index=False)
```

## Quality Metrics

### Quality Score Calculation

```python
peak_score = num_matched_peaks / min_explained_peaks
intensity_score = matched_intensity_pct / 100
quality_score = 0.3 * peak_score + 0.7 * intensity_score
```

**Classification:**
- **good**: quality_score ≥ 0.7 AND num_matched_peaks ≥ min_peaks
- **uncertain**: 0.4 ≤ quality_score < 0.7
- **bad**: quality_score < 0.4

### Cosine Similarity

Modified cosine score with sqrt intensity scaling:

1. Normalize intensities: `intensity_norm = intensity^0.5 / ||intensity^0.5||`
2. Match peaks within tolerance
3. Calculate: `cosine_score = Σ(intensity1[i] × intensity2[i])` for matched peaks

**Interpretation:**
- **0.9-1.0**: Nearly identical spectra
- **0.7-0.9**: Very similar (likely isomers)
- **0.5-0.7**: Moderately similar
- **<0.5**: Dissimilar

## Troubleshooting

### No valid spectra found

**Cause:** MGF missing SMILES or FORMULA fields

**Solution:**
```bash
# Check MGF has required fields
grep -E "SMILES|FORMULA" your_file.mgf | head -20
```

### Cosine similarity very slow

**Cause:** Too many spectra (O(n²) complexity)

**Solutions:**
1. Increase `--min_similarity 0.7`
2. Reduce `--cosine_top_n 5`
3. Split into smaller batches
4. Use more memory: `--max_memory 32.GB`

### Low annotation quality

**Causes:**
1. Incorrect formulas in reference
2. Poor fragmentation
3. Thresholds too strict

**Solutions:**
1. Check `mass_error_ppm` in reference.tsv
2. Lower thresholds: `--quality_threshold_good 0.5`
3. Review fragment generation

## Performance

| Dataset Size | Runtime | Memory |Cosine Similarity Time |
|-------------|---------|--------|----------------------|
| 100 spectra | 5-10 min | 2-4 GB | 1-2 min |
| 1,000 spectra | 20-40 min | 4-8 GB | 10-20 min |
| 10,000 spectra | 2-4 hours | 8-16 GB | 2-3 hours |
| 100,000 spectra | 10-20 hours | 32-64 GB | 10-15 hours |

**Note:** Cosine similarity calculation dominates runtime for large datasets

## Citation

If you use this pipeline, please cite:

**GNPS:**
- Wang, M. et al. (2016). Sharing and community curation of mass spectrometry data with Global Natural Products Social Molecular Networking. *Nature Biotechnology*, 34(8), 828-837.

## Support

- **Issues**: [GitHub Issues](https://github.com/yourusername/gnps-annotation-pipeline/issues)
- **Documentation**: See `GNPS_MODE.md` for detailed documentation
- **Examples**: See `examples/` directory

## License

MIT License - see LICENSE file for details

## Authors

- **Bigy Ambat** - Project Lead
- **Contributors** - See CONTRIBUTORS.md

## Version History

- **3.0** (2025-11-12) - GNPS-only mode, removed MSBuddy dependency
- **2.1** (Previous) - MSBuddy prediction mode
