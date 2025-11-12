# GNPS Reference Library Mode - Implementation Summary

## Overview

This document summarizes the implementation of GNPS Reference Library mode for the MSBuddy-Annotation-Pipeline. The update enables the pipeline to work with pre-annotated GNPS spectral library data for quality assessment and spectral similarity analysis.

## Implementation Date

2025-11-12

## Key Changes

### 1. New Python Scripts

Four new Python scripts were added to `bin/`:

#### a. `parse_gnps_reference.py` (602 lines)
**Purpose:** Extract reference annotations from GNPS MGF files

**Key Functions:**
- `parse_gnps_mgf()`: Loads MGF and extracts SMILES, InChI, formulas
- `parse_molecular_formula()`: Parses formula strings into element counts
- `calculate_formula_mass()`: Calculates monoisotopic mass
- `generate_theoretical_fragments()`: Creates theoretical fragment library
- `validate_annotations()`: Reports annotation completeness statistics

**Input:** GNPS MGF file with SMILES/formula metadata
**Output:** TSV with parsed reference annotations

#### b. `annotate_peaks_gnps.py` (819 lines)
**Purpose:** Annotate fragment peaks using reference formulas

**Key Functions:**
- `load_mgf_spectra()`: Loads spectrum data
- `generate_fragment_library()`: Creates theoretical fragments from formula
- `match_peaks_to_fragments()`: Matches observed peaks to theoretical fragments
- `calculate_quality_score()`: Computes overall quality score (0-1)
- `classify_quality()`: Classifies as good/uncertain/bad/no_annotation

**Input:** MGF file + reference annotations TSV
**Output:** Annotated peaks TSV with quality metrics

#### c. `calculate_cosine_similarity.py` (562 lines)
**Purpose:** Calculate pairwise spectral similarities

**Key Functions:**
- `normalize_spectrum()`: Normalizes intensities for comparison
- `calculate_cosine_similarity()`: Computes modified cosine score
- `calculate_pairwise_similarities()`: All-vs-all similarity matrix
- `find_top_neighbors()`: Identifies most similar spectra
- `enhance_annotations_with_similarity()`: Adds similarity metrics to annotations

**Input:** MGF file + annotated peaks TSV
**Output:** Enhanced annotations TSV + neighbors TSV

#### d. `generate_qc_report_gnps.py` (977 lines)
**Purpose:** Generate HTML QC report for GNPS workflow

**Key Functions:**
- `plot_annotation_quality_distribution()`: Quality classification bar chart
- `plot_quality_score_distribution()`: Score histogram
- `plot_matched_peaks_distribution()`: Peak matching histogram
- `plot_similarity_distribution()`: Cosine similarity histogram
- `plot_quality_vs_similarity()`: Scatter plot of quality vs. similarity
- `generate_html_report()`: Creates self-contained HTML with embedded plots

**Input:** Enhanced annotations TSV
**Output:** Interactive HTML report

### 2. Updated Nextflow Workflow (`main.nf`)

**Changes:**
- Added version bump to 3.0
- Added mode-based header printing (MSBuddy vs. GNPS)
- Added 4 new GNPS processes:
  - `PARSE_GNPS_REFERENCE`
  - `ANNOTATE_PEAKS_GNPS`
  - `CALCULATE_COSINE_SIMILARITY`
  - `GENERATE_QC_GNPS`
- Implemented conditional workflow branching based on `params.gnps_mode`

**GNPS Workflow Chain:**
```
MGF → PARSE_GNPS_REFERENCE → ANNOTATE_PEAKS_GNPS →
CALCULATE_COSINE_SIMILARITY → GENERATE_QC_GNPS → HTML Report
```

### 3. Updated Configuration (`nextflow.config`)

**New Parameters:**
```groovy
params {
    gnps_mode               = false     // Enable GNPS mode
    gnps_mz_tol             = 0.01      // m/z tolerance (Da)
    gnps_ppm_tol            = 20.0      // PPM tolerance
    gnps_min_peaks          = 3         // Min peaks in spectrum
    gnps_cosine_top_n       = 10        // Top N neighbors
    gnps_min_similarity     = 0.5       // Min similarity score
}
```

**New Process Resources:**
```groovy
PARSE_GNPS_REFERENCE:       2 CPU, 4 GB, 1h
ANNOTATE_PEAKS_GNPS:        4 CPU, 8 GB, 2h
CALCULATE_COSINE_SIMILARITY: 4 CPU, 16 GB, 4h
GENERATE_QC_GNPS:           1 CPU, 4 GB, 1h
```

### 4. Documentation

**New Files:**
- `GNPS_MODE.md`: Comprehensive user guide (430 lines)
  - Pipeline workflow explanation
  - Input format requirements
  - Parameter descriptions
  - Output file descriptions
  - Use case examples
  - Troubleshooting guide
  - Comparison with MSBuddy mode

- `GNPS_IMPLEMENTATION_SUMMARY.md`: This file

## Architecture

### Data Flow

```
┌─────────────────────────────────────────────────────────────┐
│ GNPS MGF Input                                              │
│ - Contains: SMILES, Formula, SPECTRUMID                     │
│ - Fields: SMILES, FORMULA, PEPMASS, fragment peaks          │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│ PARSE_GNPS_REFERENCE                                        │
│ - Extracts reference annotations                            │
│ - Validates completeness                                    │
│ - Calculates formula mass                                   │
│ Output: reference.tsv                                       │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│ ANNOTATE_PEAKS_GNPS                                         │
│ - Generates theoretical fragments from formula              │
│ - Matches peaks to fragments                                │
│ - Calculates quality metrics                                │
│ - Classifies: good/uncertain/bad                            │
│ Output: annotated.tsv                                       │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│ CALCULATE_COSINE_SIMILARITY                                 │
│ - Normalizes spectra                                        │
│ - Computes pairwise similarities                            │
│ - Finds top N neighbors                                     │
│ - Adds similarity metrics                                   │
│ Output: with_similarity.tsv, neighbors.tsv                  │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│ GENERATE_QC_GNPS                                            │
│ - Generates 8 visualization plots                           │
│ - Creates summary statistics                                │
│ - Builds self-contained HTML                                │
│ Output: qc_report.html                                      │
└─────────────────────────────────────────────────────────────┘
```

### Quality Metrics

**Quality Score Calculation:**
```python
peak_score = num_matched_peaks / min_explained_peaks
intensity_score = matched_intensity_pct / 100
quality_score = 0.3 * peak_score + 0.7 * intensity_score
```

**Quality Classification:**
- **good**: quality_score ≥ 0.7 AND num_matched ≥ min_peaks
- **uncertain**: 0.4 ≤ quality_score < 0.7
- **bad**: quality_score < 0.4
- **no_annotation**: num_matched = 0

**Cosine Similarity:**
```python
# Normalize spectra
intensity_normalized = intensity^0.5 / ||intensity^0.5||

# Match peaks within tolerance
matches = find_peaks_within_tolerance(mz1, mz2, tolerance)

# Calculate dot product
cosine_score = sum(intensity1[i] * intensity2[i] for i in matches)
```

## Output Files

### Directory Structure

```
results/
├── gnps_reference/
│   └── {sample}_reference.tsv          (Parsed annotations)
├── gnps_annotations/
│   └── {sample}_annotated.tsv          (Peak annotations + quality)
├── gnps_similarity/
│   ├── {sample}_with_similarity.tsv    (Enhanced with cosine scores)
│   └── {sample}_neighbors.tsv          (Top N similar spectra)
├── qc_reports/
│   └── {sample}_qc_report.html         (Interactive report)
└── pipeline_info/
    ├── execution_timeline.html
    ├── execution_report.html
    └── execution_trace.txt
```

### Key Output Columns

**Reference TSV:**
- `spectrum_id`, `smiles`, `inchi`, `inchikey`, `formula`, `precursor_mz`
- `charge`, `adduct`, `num_peaks`, `formula_mass`, `compound_name`
- `has_structure`, `has_formula`, `annotation_complete`

**Annotated TSV:**
- All reference columns plus:
- `num_matched_peaks`, `matched_peaks_pct`, `matched_intensity_pct`
- `quality_score`, `annotation_quality`, `num_theoretical_fragments`

**Similarity TSV:**
- All annotated columns plus:
- `num_similar_neighbors`, `max_similarity_score`, `avg_similarity_score`
- `top_neighbor`, `top_neighbor_score`

**Neighbors TSV:**
- `query_spectrum`, `match_spectrum`, `rank`, `cosine_score`
- `num_matching_peaks`, `precursor_mz_diff`

## Usage Examples

### Basic GNPS Mode

```bash
nextflow run main.nf \
    --input 'gnps_data/*.mgf' \
    --outdir results_gnps \
    --gnps_mode true \
    -profile docker
```

### High-Resolution Data

```bash
nextflow run main.nf \
    --input 'orbitrap_data/*.mgf' \
    --outdir results \
    --gnps_mode true \
    --gnps_mz_tol 0.005 \
    --gnps_ppm_tol 5.0 \
    -profile docker
```

### Large Dataset (HPC)

```bash
nextflow run main.nf \
    --input 'large_library/*.mgf' \
    --outdir results \
    --gnps_mode true \
    --gnps_min_similarity 0.7 \
    --max_memory 32.GB \
    -profile slurm,singularity
```

## Testing Checklist

- [ ] Test with sample GNPS MGF file
- [ ] Verify reference parsing extracts SMILES correctly
- [ ] Check peak annotation generates quality scores
- [ ] Confirm cosine similarity calculation completes
- [ ] Validate HTML report generation
- [ ] Test with missing SMILES/formula (graceful failure)
- [ ] Test with large dataset (>1000 spectra)
- [ ] Verify Docker container compatibility
- [ ] Check Singularity compatibility
- [ ] Test on HPC cluster (SLURM)

## Performance Considerations

### Computational Complexity

- **Reference Parsing:** O(n) - linear in number of spectra
- **Peak Annotation:** O(n × m) - n spectra, m theoretical fragments
- **Cosine Similarity:** O(n² × p) - n² pairwise comparisons, p peaks
  - **Bottleneck for large datasets!**
  - Optimization: Use `--gnps_min_similarity` threshold
  - Optimization: Reduce `--gnps_cosine_top_n`

### Memory Requirements

| Process | Typical Memory | Large Dataset |
|---------|---------------|---------------|
| PARSE_GNPS_REFERENCE | 2-4 GB | 4-8 GB |
| ANNOTATE_PEAKS_GNPS | 4-8 GB | 8-16 GB |
| CALCULATE_COSINE_SIMILARITY | 8-16 GB | 32-64 GB |
| GENERATE_QC_GNPS | 2-4 GB | 4-8 GB |

### Estimated Runtime

| Dataset Size | Total Runtime | Cosine Similarity |
|-------------|--------------|-------------------|
| 100 spectra | 5-10 min | 1-2 min |
| 1,000 spectra | 20-40 min | 10-20 min |
| 10,000 spectra | 2-4 hours | 2-3 hours |
| 100,000 spectra | 10-20 hours | 10-15 hours |

**Note:** Cosine similarity dominates runtime for large datasets (O(n²))

## Future Enhancements

### Planned Features

1. **SQLite Database Integration**
   - Store precomputed similarities
   - Enable fast querying
   - Support incremental updates

2. **Metadata Enrichment**
   - Import collision energy, CCS, RT from external DB
   - Impute missing metadata from similar spectra
   - Weighted average based on cosine scores

3. **Spectral Quality Reranking**
   - Combine peak quality + cosine similarity
   - ML-based quality prediction
   - Integration with GNPS workflows

4. **Parallel Cosine Calculation**
   - GPU acceleration for large datasets
   - Chunked processing for memory efficiency
   - Distributed computing support

5. **Advanced Fragmentation Models**
   - Structure-based fragment prediction (CFM-ID style)
   - Substructure annotation (MAGMa style)
   - Neutral loss libraries

### Known Limitations

1. **Theoretical Fragmentation:**
   - Current implementation uses simple neutral losses
   - Doesn't account for structure-specific fragmentation
   - No bond dissociation energy considerations

2. **Cosine Similarity:**
   - Computationally expensive for large datasets
   - Requires all-vs-all comparison (no indexing yet)
   - Memory-intensive for >10,000 spectra

3. **Quality Thresholds:**
   - Fixed thresholds may not work for all compound classes
   - No automatic threshold optimization
   - Doesn't account for instrument-specific biases

## Integration with Existing Workflow

### Backward Compatibility

- **MSBuddy mode still default:** `gnps_mode = false`
- **All existing functionality preserved**
- **Separate output directories** prevent conflicts
- **Independent process definitions** for each mode

### Code Reuse

**Shared Functions:**
- MGF parsing logic (from `analyze_peak_explanation.py`)
- Spectrum ID extraction (priority order: SPECTRUMID → TITLE → index)
- HTML report styling (adapted from `generate_qc_report.py`)

**Mode-Specific:**
- Formula prediction (MSBuddy) vs. reference formulas (GNPS)
- FDR-based quality (MSBuddy) vs. peak-based quality (GNPS)
- No similarity calculation (MSBuddy) vs. cosine similarity (GNPS)

## Dependencies

### New Python Libraries

All dependencies already present in `environment.yml`:
- `pyteomics >= 4.5.0` (MGF parsing)
- `pandas >= 1.5.0` (Data manipulation)
- `numpy >= 1.23.0` (Numerical operations)
- `matplotlib >= 3.6.0` (Plotting)
- `seaborn >= 0.12.0` (Statistical plots)

**Note:** No additional dependencies needed!

## Version Control

### Modified Files

```
main.nf                           (v2.0 → v3.0, +100 lines)
nextflow.config                   (v2.0 → v3.0, +35 lines)
```

### New Files

```
bin/parse_gnps_reference.py       (602 lines)
bin/annotate_peaks_gnps.py        (819 lines)
bin/calculate_cosine_similarity.py (562 lines)
bin/generate_qc_report_gnps.py    (977 lines)
GNPS_MODE.md                      (430 lines)
GNPS_IMPLEMENTATION_SUMMARY.md    (this file)
```

**Total New Code:** ~3,400 lines

## Contributors

- Implementation: Claude (Anthropic)
- Project Lead: Bigy Ambat
- Date: 2025-11-12

## License

MIT License (same as parent project)

---

## Quick Reference

### Enable GNPS Mode

```bash
--gnps_mode true
```

### Key Parameters

```bash
--gnps_mz_tol 0.01              # m/z tolerance (Da)
--gnps_ppm_tol 20.0             # PPM tolerance
--gnps_min_similarity 0.5       # Min cosine score
--gnps_cosine_top_n 10          # Top N neighbors
```

### Output Locations

```bash
results/gnps_annotations/       # Peak annotations
results/gnps_similarity/        # Cosine scores
results/qc_reports/             # HTML reports
```

### Quality Filter

```python
df = pd.read_csv('results/gnps_similarity/data_with_similarity.tsv', sep='\t')
good = df[df['annotation_quality'] == 'good']
```
