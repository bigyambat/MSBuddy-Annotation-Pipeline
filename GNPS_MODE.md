# GNPS Reference Library Mode

## Overview

The GNPS Reference Library mode enables annotation quality assessment for GNPS spectral libraries. Unlike the standard MSBuddy prediction mode, this mode uses **reference SMILES and molecular formulas** already present in GNPS MGF files to:

1. **Annotate explained peaks** using reference molecular formulas
2. **Separate good & bad annotations** based on peak explanation quality
3. **Calculate spectral similarity** within the dataset using cosine scores
4. **Generate quality metrics** for spectral library refinement

This mode is designed for researchers who want to:
- Assess quality of GNPS reference library annotations
- Identify high-confidence vs. low-confidence spectral annotations
- Find spectrally similar compounds for metadata imputation
- Prepare clean training data for machine learning models

## Pipeline Workflow (GNPS Mode)

```
GNPS MGF Files (with SMILES & Formulas)
    ↓
[PARSE_GNPS_REFERENCE]
    ↓ Extract reference annotations
    ↓
Reference Annotations TSV
    ↓
[ANNOTATE_PEAKS_GNPS]
    ↓ Match theoretical fragments to observed peaks
    ↓
Annotated Peaks TSV (with quality metrics)
    ↓
[CALCULATE_COSINE_SIMILARITY]
    ↓ Compute pairwise spectral similarities
    ↓
Enhanced Annotations TSV (with similarity data)
    ↓
[GENERATE_QC_GNPS]
    ↓
HTML QC Reports
```

## Input Requirements

### GNPS MGF Format

Your MGF files should contain GNPS reference library data with the following fields:

```
BEGIN IONS
TITLE=Spectrum_Name
SPECTRUMID=CCMSLIB00000001234        ← GNPS library ID
SMILES=CC(C)CC1=CC=C(C=C1)C(C)C      ← Required: SMILES structure
FORMULA=C15H24                        ← Required: Molecular formula
INCHI=InChI=1S/C15H24/...            ← Optional: InChI
INCHIKEY=XMGQYMWWDOXHJM-UHFFFAOYSA-N ← Optional: InChIKey
NAME=Compound Name                    ← Optional: Compound name
PEPMASS=205.1956                      ← Precursor m/z
CHARGE=1+                             ← Charge state
123.4567 1000.0                       ← Fragment peaks
234.5678 500.0
END IONS
```

**Required Fields:**
- `SMILES`: Chemical structure in SMILES format
- `FORMULA` or `MOLECULARFORMULA`: Molecular formula
- `PEPMASS`: Precursor m/z value
- Fragment peak list (m/z intensity pairs)

**Optional Fields:**
- `SPECTRUMID`: GNPS library identifier
- `INCHI`/`INCHIKEY`: Chemical identifiers
- `NAME` or `COMPOUND_NAME`: Compound name
- `COLLISION_ENERGY`: Collision energy used
- `RTINSECONDS`: Retention time in seconds

## Quick Start

### Enable GNPS Mode

To run the pipeline in GNPS mode, add the `--gnps_mode true` flag:

```bash
nextflow run main.nf \
    --input 'gnps_data/*.mgf' \
    --outdir results_gnps \
    --gnps_mode true \
    -profile docker
```

### Complete Example

```bash
nextflow run main.nf \
    --input 'GNPS_library/*.mgf' \
    --outdir gnps_results \
    --gnps_mode true \
    --gnps_mz_tol 0.01 \
    --gnps_ppm_tol 20 \
    --gnps_min_similarity 0.6 \
    --gnps_cosine_top_n 10 \
    -profile docker
```

## Parameters

### GNPS-Specific Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--gnps_mode` | `false` | Enable GNPS reference library mode |
| `--gnps_mz_tol` | `0.01` | m/z tolerance for peak matching (Da) |
| `--gnps_ppm_tol` | `20.0` | PPM tolerance for peak matching |
| `--gnps_min_peaks` | `3` | Minimum peaks required in spectrum |
| `--gnps_cosine_top_n` | `10` | Number of top similar neighbors to report |
| `--gnps_min_similarity` | `0.5` | Minimum cosine similarity score to report |

### Quality Threshold Parameters (Shared)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--score_threshold` | `0.7` | Minimum quality score for "good" classification |
| `--explained_peaks_pct_threshold` | `0.5` | Minimum explained peaks percentage |
| `--explained_intensity_pct_threshold` | `0.6` | Minimum explained intensity percentage |
| `--min_explained_peaks` | `3` | Minimum number of explained peaks |

## Output Structure

### GNPS Mode Output Directory

```
results_gnps/
├── gnps_reference/                       ← Parsed reference annotations
│   ├── sample1_reference.tsv
│   └── sample2_reference.tsv
├── gnps_annotations/                     ← Peak annotations with quality metrics
│   ├── sample1_annotated.tsv
│   └── sample2_annotated.tsv
├── gnps_similarity/                      ← Cosine similarity results
│   ├── sample1_with_similarity.tsv
│   ├── sample1_neighbors.tsv
│   ├── sample2_with_similarity.tsv
│   └── sample2_neighbors.tsv
├── qc_reports/                           ← HTML QC reports
│   ├── sample1_qc_report.html
│   └── sample2_qc_report.html
└── pipeline_info/                        ← Nextflow execution logs
    ├── execution_timeline.html
    ├── execution_report.html
    └── execution_trace.txt
```

### Output File Descriptions

#### 1. Reference Annotations TSV

**File:** `gnps_reference/{sample}_reference.tsv`

Contains parsed reference data from GNPS MGF:

| Column | Description |
|--------|-------------|
| `spectrum_id` | Spectrum identifier (SPECTRUMID or TITLE) |
| `smiles` | SMILES structure |
| `inchi` | InChI identifier |
| `inchikey` | InChIKey |
| `formula` | Molecular formula |
| `precursor_mz` | Precursor m/z |
| `charge` | Ion charge state |
| `adduct` | Adduct type |
| `num_peaks` | Number of peaks in spectrum |
| `formula_mass` | Calculated mass from formula |
| `compound_name` | Compound name (if available) |
| `library_id` | GNPS library ID |
| `has_structure` | Boolean: has SMILES or InChI |
| `has_formula` | Boolean: has molecular formula |
| `annotation_complete` | Boolean: has both SMILES and formula |

#### 2. Annotated Peaks TSV

**File:** `gnps_annotations/{sample}_annotated.tsv`

Contains peak annotation results with quality metrics:

| Column | Description |
|--------|-------------|
| `spectrum_id` | Spectrum identifier |
| `smiles` | SMILES structure |
| `formula` | Molecular formula |
| `compound_name` | Compound name |
| `precursor_mz` | Precursor m/z |
| `num_peaks` | Total peaks in spectrum |
| `num_matched_peaks` | Number of peaks matched to theoretical fragments |
| `matched_peaks_pct` | Percentage of peaks matched (0-100) |
| `matched_intensity_pct` | Percentage of intensity matched (0-100) |
| `quality_score` | Overall quality score (0-1) |
| `annotation_quality` | Classification: good/uncertain/bad/no_annotation |
| `num_theoretical_fragments` | Number of theoretical fragments generated |

**Quality Classifications:**
- **good**: `quality_score >= 0.7` AND `num_matched_peaks >= min_peaks`
- **uncertain**: `0.4 <= quality_score < 0.7`
- **bad**: `quality_score < 0.4`
- **no_annotation**: No matched peaks

#### 3. Similarity-Enhanced Annotations TSV

**File:** `gnps_similarity/{sample}_with_similarity.tsv`

Enhanced annotations with spectral similarity metrics (all columns from annotated TSV plus):

| Column | Description |
|--------|-------------|
| `num_similar_neighbors` | Number of spectrally similar neighbors found |
| `max_similarity_score` | Maximum cosine similarity score (0-1) |
| `avg_similarity_score` | Average similarity score across neighbors |
| `top_neighbor` | Spectrum ID of most similar neighbor |
| `top_neighbor_score` | Cosine score of top neighbor |

#### 4. Top Neighbors TSV

**File:** `gnps_similarity/{sample}_neighbors.tsv`

Detailed similarity relationships:

| Column | Description |
|--------|-------------|
| `query_spectrum` | Query spectrum ID |
| `match_spectrum` | Matched spectrum ID |
| `rank` | Similarity rank (1 = most similar) |
| `cosine_score` | Cosine similarity score (0-1) |
| `num_matching_peaks` | Number of peaks that matched |
| `precursor_mz_diff` | Difference in precursor m/z |

#### 5. QC Report HTML

**File:** `qc_reports/{sample}_qc_report.html`

Interactive HTML report with visualizations:

**Included Plots:**
1. **Annotation Quality Distribution** - Bar chart of good/uncertain/bad/no_annotation counts
2. **Quality Score Distribution** - Histogram of quality scores with thresholds
3. **Matched Peaks Distribution** - Histogram of matched peak percentages
4. **Matched Intensity Distribution** - Histogram of matched intensity percentages
5. **Spectral Similarity Distribution** - Histogram of max cosine similarity scores
6. **Quality vs. Similarity Scatter Plot** - Correlation between peak quality and spectral similarity
7. **Similar Neighbors Distribution** - Histogram of neighbor counts per spectrum
8. **Adduct Frequency** - Top 20 most common adduct types

## Understanding Quality Metrics

### Quality Score Calculation

The quality score (0-1) combines two components:

```python
peak_score = num_matched_peaks / min_explained_peaks
intensity_score = matched_intensity_pct / 100

quality_score = 0.3 * peak_score + 0.7 * intensity_score
```

**Interpretation:**
- **0.7-1.0**: Excellent annotation quality (good)
- **0.4-0.7**: Moderate quality (uncertain)
- **0.0-0.4**: Poor quality (bad)

### Cosine Similarity

Spectral similarity is calculated using the modified cosine score:

1. **Normalize intensities**: Apply square-root scaling and L2 normalization
2. **Match peaks**: Find peaks within m/z tolerance
3. **Calculate dot product**: Sum of (intensity1 × intensity2) for matched peaks

**Interpretation:**
- **0.9-1.0**: Nearly identical spectra
- **0.7-0.9**: Very similar spectra (likely isomers or close analogs)
- **0.5-0.7**: Moderately similar (shared structural features)
- **0.0-0.5**: Dissimilar spectra

## Use Cases

### 1. Library Quality Assessment

Identify high-quality vs. low-quality reference annotations:

```bash
nextflow run main.nf \
    --input 'gnps_library/*.mgf' \
    --outdir library_qc \
    --gnps_mode true \
    --score_threshold 0.8 \
    -profile docker
```

**Filter good annotations:**
```python
import pandas as pd

df = pd.read_csv('library_qc/gnps_similarity/library_with_similarity.tsv', sep='\t')
good_annotations = df[df['annotation_quality'] == 'good']
good_annotations.to_csv('high_quality_library.tsv', sep='\t', index=False)
```

### 2. Spectral Clustering

Find spectral clusters for molecular network analysis:

```bash
nextflow run main.nf \
    --input 'gnps_data/*.mgf' \
    --outdir clustering \
    --gnps_mode true \
    --gnps_min_similarity 0.7 \
    --gnps_cosine_top_n 20 \
    -profile docker
```

**Analyze clusters:**
```python
import pandas as pd
import networkx as nx

neighbors = pd.read_csv('clustering/gnps_similarity/data_neighbors.tsv', sep='\t')

# Build network
G = nx.from_pandas_edgelist(
    neighbors,
    source='query_spectrum',
    target='match_spectrum',
    edge_attr='cosine_score'
)

# Find connected components (clusters)
clusters = list(nx.connected_components(G))
print(f"Found {len(clusters)} spectral clusters")
```

### 3. Training Data Preparation for ML

Extract clean, high-quality annotations for machine learning:

```bash
nextflow run main.nf \
    --input 'training_data/*.mgf' \
    --outdir ml_prep \
    --gnps_mode true \
    --score_threshold 0.8 \
    --min_explained_peaks 5 \
    -profile docker
```

**Filter for ML training:**
```python
import pandas as pd

df = pd.read_csv('ml_prep/gnps_similarity/data_with_similarity.tsv', sep='\t')

# High-quality annotations only
ml_data = df[
    (df['annotation_quality'] == 'good') &
    (df['matched_intensity_pct'] > 70) &
    (df['num_matched_peaks'] >= 5) &
    (df['smiles'].notna())
]

ml_data[['spectrum_id', 'smiles', 'formula', 'quality_score']].to_csv(
    'ml_training_set.tsv', sep='\t', index=False
)
```

### 4. Metadata Imputation from Neighbors

Use similar spectra to impute missing metadata:

```python
import pandas as pd
import numpy as np

annotations = pd.read_csv('results/gnps_similarity/data_with_similarity.tsv', sep='\t')
neighbors = pd.read_csv('results/gnps_similarity/data_neighbors.tsv', sep='\t')

# Find spectra missing collision energy
missing_ce = annotations[annotations['collision_energy'].isna()]

# Impute from similar neighbors
for idx, row in missing_ce.iterrows():
    spec_id = row['spectrum_id']

    # Get neighbors
    similar = neighbors[neighbors['query_spectrum'] == spec_id]

    # Find neighbors with known CE
    neighbor_ids = similar['match_spectrum'].tolist()
    neighbor_data = annotations[
        (annotations['spectrum_id'].isin(neighbor_ids)) &
        (annotations['collision_energy'].notna())
    ]

    if len(neighbor_data) > 0:
        # Weight by similarity score
        weights = similar[similar['match_spectrum'].isin(neighbor_data['spectrum_id'])]['cosine_score']
        imputed_ce = np.average(neighbor_data['collision_energy'], weights=weights)
        annotations.at[idx, 'collision_energy'] = imputed_ce
```

## Advanced Configuration

### High-Resolution MS Data

For high-resolution instruments (Orbitrap, FTICR):

```bash
nextflow run main.nf \
    --gnps_mode true \
    --gnps_mz_tol 0.005 \
    --gnps_ppm_tol 5.0 \
    -profile docker
```

### Low-Resolution MS Data

For lower-resolution instruments (ion trap, TOF):

```bash
nextflow run main.nf \
    --gnps_mode true \
    --gnps_mz_tol 0.05 \
    --gnps_ppm_tol 50.0 \
    -profile docker
```

### Large Dataset Processing

For very large datasets (>10,000 spectra):

```bash
nextflow run main.nf \
    --gnps_mode true \
    --gnps_min_similarity 0.7 \
    --gnps_cosine_top_n 5 \
    --max_memory 32.GB \
    -profile slurm
```

## Troubleshooting

### Issue: "No valid spectra found"

**Cause:** MGF file missing required SMILES or FORMULA fields

**Solution:** Check your MGF file has required fields:
```bash
grep -E "SMILES|FORMULA" your_file.mgf | head -20
```

### Issue: "Cosine similarity calculation very slow"

**Cause:** Too many spectra for all-vs-all comparison

**Solutions:**
1. Increase `--gnps_min_similarity` threshold (e.g., 0.7)
2. Reduce `--gnps_cosine_top_n` (e.g., 5)
3. Process in smaller batches
4. Use more memory (`--max_memory 32.GB`)

### Issue: "Low annotation quality across board"

**Possible Causes:**
1. Incorrect molecular formulas in reference
2. Poor fragmentation patterns
3. Thresholds too strict

**Solutions:**
1. Check `mass_error_ppm` in reference TSV for formula accuracy
2. Lower quality thresholds: `--score_threshold 0.5`
3. Review theoretical fragment generation in `parse_gnps_reference.py`

## Comparison: GNPS Mode vs. MSBuddy Mode

| Feature | MSBuddy Mode | GNPS Mode |
|---------|-------------|-----------|
| **Input** | Raw MGF (no formulas needed) | GNPS MGF (requires SMILES/formula) |
| **Formula Source** | Predicted by MSBuddy | From reference library |
| **Primary Quality Metric** | FDR (False Discovery Rate) | Peak explanation quality |
| **Fragment Annotation** | MSBuddy ML models | Theoretical fragmentation |
| **Spectral Similarity** | Not calculated | Cosine similarity computed |
| **Best For** | Unknown compound identification | Library curation & validation |
| **Speed** | Slower (ML inference) | Faster (no prediction) |
| **Accuracy** | High (95% at FDR < 0.05) | Depends on reference quality |

## Citation

If you use this pipeline in your research, please cite:

**For MSBuddy:**
- Original MSBuddy paper citation here

**For GNPS:**
- Wang, M. et al. (2016). Sharing and community curation of mass spectrometry data with Global Natural Products Social Molecular Networking. Nature Biotechnology, 34(8), 828-837.

## Support

For issues, questions, or feature requests, please open an issue on GitHub or contact the maintainers.

## License

MIT License - see LICENSE file for details
