#!/usr/bin/env python3
"""
Analyze peak explanation from MSBuddy annotations and classify annotations as good/bad.

This script:
1. Parses MSBuddy annotation results
2. Loads corresponding MGF spectra
3. Calculates peak explanation statistics for each spectrum
4. Classifies annotations as good/bad based on configurable thresholds
5. Outputs enhanced TSV with all quality metrics
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
from pyteomics import mgf
import numpy as np
from typing import Dict, List, Tuple, Optional


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Analyze peak explanation and classify MSBuddy annotations'
    )
    parser.add_argument(
        '--mgf',
        type=str,
        required=True,
        help='Path to input MGF file'
    )
    parser.add_argument(
        '--msbuddy_tsv',
        type=str,
        required=True,
        help='Path to MSBuddy annotation TSV file'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Path to output enhanced TSV file'
    )
    parser.add_argument(
        '--score_threshold',
        type=float,
        default=0.7,
        help='Minimum MSBuddy score for good annotation (default: 0.7)'
    )
    parser.add_argument(
        '--explained_peaks_pct_threshold',
        type=float,
        default=0.5,
        help='Minimum percentage of explained peaks for good annotation (default: 0.5)'
    )
    parser.add_argument(
        '--explained_intensity_pct_threshold',
        type=float,
        default=0.6,
        help='Minimum percentage of explained intensity for good annotation (default: 0.6)'
    )
    parser.add_argument(
        '--mass_error_threshold',
        type=float,
        default=5.0,
        help='Maximum mass error (ppm) for good annotation (default: 5.0)'
    )
    parser.add_argument(
        '--min_explained_peaks',
        type=int,
        default=3,
        help='Minimum number of explained peaks for good annotation (default: 3)'
    )

    return parser.parse_args()


def load_mgf_spectra(mgf_path: str) -> Dict:
    """
    Load spectra from MGF file into a dictionary keyed by spectrum title/ID.

    Args:
        mgf_path: Path to MGF file

    Returns:
        Dictionary mapping spectrum ID to spectrum data
    """
    spectra_dict = {}

    try:
        with mgf.read(mgf_path) as reader:
            for idx, spectrum in enumerate(reader):
                # Skip None or empty spectra
                if spectrum is None:
                    continue

                # Skip spectra without required fields
                if 'params' not in spectrum or 'm/z array' not in spectrum or 'intensity array' not in spectrum:
                    continue

                # Try to get spectrum ID from various possible fields
                spec_id = spectrum.get('params', {}).get('title',
                          spectrum.get('params', {}).get('scans',
                          spectrum.get('params', {}).get('spectrum_id', f'spectrum_{idx}')))

                # Handle case where title might be directly in spectrum dict
                if spec_id == f'spectrum_{idx}' and 'title' in spectrum:
                    spec_id = spectrum['title']

                spectra_dict[spec_id] = {
                    'mz': spectrum['m/z array'],
                    'intensity': spectrum['intensity array'],
                    'precursor_mz': spectrum['params'].get('pepmass', [0])[0] if isinstance(spectrum['params'].get('pepmass'), list) else spectrum['params'].get('pepmass', 0),
                    'charge': spectrum['params'].get('charge', [1])[0] if isinstance(spectrum['params'].get('charge'), list) else spectrum['params'].get('charge', 1),
                }

    except Exception as e:
        print(f"Error loading MGF file: {e}", file=sys.stderr)
        raise

    return spectra_dict


def load_msbuddy_annotations(tsv_path: str) -> pd.DataFrame:
    """
    Load MSBuddy annotation results from TSV file.

    Args:
        tsv_path: Path to MSBuddy TSV file

    Returns:
        DataFrame with annotation results
    """
    try:
        df = pd.read_csv(tsv_path, sep='\t')
        return df
    except Exception as e:
        print(f"Error loading MSBuddy TSV file: {e}", file=sys.stderr)
        raise


def parse_fragment_annotations(annotation_str: Optional[str]) -> List[Dict]:
    """
    Parse MSBuddy fragment annotation string to extract explained peaks.

    MSBuddy typically outputs fragment annotations in formats like:
    - "123.456:formula1;234.567:formula2"
    - JSON format
    - or other structured formats

    Args:
        annotation_str: Fragment annotation string from MSBuddy

    Returns:
        List of dictionaries with fragment information
    """
    fragments = []

    if pd.isna(annotation_str) or annotation_str == '' or annotation_str is None:
        return fragments

    # Try to parse different possible formats
    # Format 1: semicolon-separated mz:formula pairs
    if ';' in str(annotation_str) and ':' in str(annotation_str):
        parts = str(annotation_str).split(';')
        for part in parts:
            if ':' in part:
                try:
                    mz_str, formula = part.split(':', 1)
                    fragments.append({
                        'mz': float(mz_str.strip()),
                        'formula': formula.strip()
                    })
                except ValueError:
                    continue

    # Format 2: comma-separated mz values (simpler format)
    elif ',' in str(annotation_str):
        parts = str(annotation_str).split(',')
        for part in parts:
            try:
                mz = float(part.strip())
                fragments.append({
                    'mz': mz,
                    'formula': 'unknown'
                })
            except ValueError:
                continue

    return fragments


def match_peaks_to_annotations(spectrum_mz: np.ndarray,
                                spectrum_intensity: np.ndarray,
                                fragment_annotations: List[Dict],
                                mz_tolerance: float = 0.01) -> Tuple[int, float, int]:
    """
    Match annotated fragments to actual peaks in the spectrum.

    Args:
        spectrum_mz: Array of m/z values
        spectrum_intensity: Array of intensity values
        fragment_annotations: List of annotated fragments
        mz_tolerance: m/z tolerance for matching (Da)

    Returns:
        Tuple of (num_explained_peaks, explained_intensity_pct, total_peaks)
    """
    if len(fragment_annotations) == 0:
        return 0, 0.0, len(spectrum_mz)

    total_peaks = len(spectrum_mz)
    total_intensity = np.sum(spectrum_intensity)

    if total_peaks == 0 or total_intensity == 0:
        return 0, 0.0, 0

    # Track which peaks are explained
    explained_peaks = set()
    explained_intensity = 0.0

    for fragment in fragment_annotations:
        frag_mz = fragment['mz']

        # Find peaks within tolerance
        mz_diffs = np.abs(spectrum_mz - frag_mz)
        matching_indices = np.where(mz_diffs <= mz_tolerance)[0]

        if len(matching_indices) > 0:
            # Take the closest match
            closest_idx = matching_indices[np.argmin(mz_diffs[matching_indices])]
            explained_peaks.add(closest_idx)
            explained_intensity += spectrum_intensity[closest_idx]

    num_explained = len(explained_peaks)
    explained_intensity_pct = (explained_intensity / total_intensity * 100) if total_intensity > 0 else 0.0

    return num_explained, explained_intensity_pct, total_peaks


def calculate_peak_statistics(row: pd.Series,
                              spectra_dict: Dict,
                              ms2_tolerance: float = 0.01) -> Dict:
    """
    Calculate peak explanation statistics for a single spectrum annotation.

    Args:
        row: DataFrame row with MSBuddy annotation
        spectra_dict: Dictionary of spectrum data
        ms2_tolerance: MS2 m/z tolerance for peak matching

    Returns:
        Dictionary with calculated statistics
    """
    stats = {
        'num_explained_peaks': 0,
        'total_peaks': 0,
        'explained_peaks_pct': 0.0,
        'explained_intensity_pct': 0.0,
        'msbuddy_score': None,
        'mass_error_ppm': None,
    }

    # Get spectrum ID (try various column names)
    spec_id = None
    for col in ['spectrum_id', 'title', 'scans', 'id', 'identifier', 'spectrum_title']:
        if col in row and pd.notna(row[col]):
            spec_id = row[col]
            break

    if spec_id is None or spec_id not in spectra_dict:
        return stats

    spectrum = spectra_dict[spec_id]
    stats['total_peaks'] = len(spectrum['mz'])

    # Extract MSBuddy score (try various column names)
    for col in ['score', 'msbuddy_score', 'confidence', 'final_score']:
        if col in row and pd.notna(row[col]):
            stats['msbuddy_score'] = float(row[col])
            break

    # Extract mass error (try various column names)
    for col in ['mass_error_ppm', 'ppm_error', 'mass_error', 'delta_ppm', 'precursor_ppm']:
        if col in row and pd.notna(row[col]):
            stats['mass_error_ppm'] = float(row[col])
            break

    # Parse fragment annotations (try various column names)
    fragment_annotations = []
    for col in ['fragments', 'fragment_annotations', 'explained_peaks', 'matched_fragments']:
        if col in row and pd.notna(row[col]):
            fragment_annotations = parse_fragment_annotations(row[col])
            if fragment_annotations:
                break

    # Calculate peak matching statistics
    if fragment_annotations:
        num_explained, explained_intensity_pct, total = match_peaks_to_annotations(
            spectrum['mz'],
            spectrum['intensity'],
            fragment_annotations,
            ms2_tolerance
        )
        stats['num_explained_peaks'] = num_explained
        stats['explained_intensity_pct'] = explained_intensity_pct
        stats['explained_peaks_pct'] = (num_explained / total * 100) if total > 0 else 0.0

    return stats


def classify_annotation(stats: Dict, thresholds: Dict) -> str:
    """
    Classify annotation as good, bad, or uncertain based on thresholds.

    Args:
        stats: Dictionary with calculated statistics
        thresholds: Dictionary with threshold values

    Returns:
        Classification: 'good', 'bad', or 'uncertain'
    """
    # Check if annotation exists (has a formula)
    if stats.get('msbuddy_score') is None:
        return 'no_annotation'

    # Count how many criteria are met
    criteria_met = 0
    total_criteria = 0

    # Criterion 1: MSBuddy score
    if stats.get('msbuddy_score') is not None:
        total_criteria += 1
        if stats['msbuddy_score'] >= thresholds['score_threshold']:
            criteria_met += 1

    # Criterion 2: Explained peaks percentage
    if stats.get('explained_peaks_pct') is not None:
        total_criteria += 1
        if stats['explained_peaks_pct'] >= thresholds['explained_peaks_pct_threshold'] * 100:
            criteria_met += 1

    # Criterion 3: Explained intensity percentage
    if stats.get('explained_intensity_pct') is not None:
        total_criteria += 1
        if stats['explained_intensity_pct'] >= thresholds['explained_intensity_pct_threshold'] * 100:
            criteria_met += 1

    # Criterion 4: Mass error
    if stats.get('mass_error_ppm') is not None:
        total_criteria += 1
        if abs(stats['mass_error_ppm']) <= thresholds['mass_error_threshold']:
            criteria_met += 1

    # Criterion 5: Minimum explained peaks count
    if stats.get('num_explained_peaks') is not None:
        total_criteria += 1
        if stats['num_explained_peaks'] >= thresholds['min_explained_peaks']:
            criteria_met += 1

    # Classification logic
    if total_criteria == 0:
        return 'uncertain'

    criteria_ratio = criteria_met / total_criteria

    if criteria_ratio >= 0.8:  # At least 80% of criteria met
        return 'good'
    elif criteria_ratio <= 0.4:  # Less than 40% of criteria met
        return 'bad'
    else:
        return 'uncertain'


def main():
    """Main execution function."""
    args = parse_args()

    print(f"Loading MGF spectra from: {args.mgf}")
    spectra_dict = load_mgf_spectra(args.mgf)
    print(f"Loaded {len(spectra_dict)} spectra")

    print(f"Loading MSBuddy annotations from: {args.msbuddy_tsv}")
    annotations_df = load_msbuddy_annotations(args.msbuddy_tsv)
    print(f"Loaded {len(annotations_df)} annotations")

    # Prepare thresholds dictionary
    thresholds = {
        'score_threshold': args.score_threshold,
        'explained_peaks_pct_threshold': args.explained_peaks_pct_threshold,
        'explained_intensity_pct_threshold': args.explained_intensity_pct_threshold,
        'mass_error_threshold': args.mass_error_threshold,
        'min_explained_peaks': args.min_explained_peaks,
    }

    # Calculate statistics for each annotation
    print("Calculating peak explanation statistics...")
    stats_list = []
    classifications = []

    for idx, row in annotations_df.iterrows():
        stats = calculate_peak_statistics(row, spectra_dict, ms2_tolerance=0.01)
        classification = classify_annotation(stats, thresholds)

        stats_list.append(stats)
        classifications.append(classification)

        if (idx + 1) % 100 == 0:
            print(f"Processed {idx + 1}/{len(annotations_df)} annotations...")

    # Add new columns to DataFrame
    annotations_df['annotation_quality'] = classifications
    annotations_df['num_explained_peaks'] = [s['num_explained_peaks'] for s in stats_list]
    annotations_df['total_peaks'] = [s['total_peaks'] for s in stats_list]
    annotations_df['explained_peaks_pct'] = [s['explained_peaks_pct'] for s in stats_list]
    annotations_df['explained_intensity_pct'] = [s['explained_intensity_pct'] for s in stats_list]

    # Add score and mass error if not already present
    if 'msbuddy_score' not in annotations_df.columns:
        annotations_df['msbuddy_score'] = [s['msbuddy_score'] for s in stats_list]
    if 'mass_error_ppm' not in annotations_df.columns:
        annotations_df['mass_error_ppm'] = [s['mass_error_ppm'] for s in stats_list]

    # Save enhanced annotations
    print(f"Saving enhanced annotations to: {args.output}")
    annotations_df.to_csv(args.output, sep='\t', index=False)

    # Print summary statistics
    print("\n=== Classification Summary ===")
    classification_counts = pd.Series(classifications).value_counts()
    for classification, count in classification_counts.items():
        pct = count / len(classifications) * 100
        print(f"{classification}: {count} ({pct:.1f}%)")

    print("\n=== Quality Statistics ===")
    if stats_list:
        explained_peaks_pcts = [s['explained_peaks_pct'] for s in stats_list if s['explained_peaks_pct'] > 0]
        if explained_peaks_pcts:
            print(f"Mean explained peaks: {np.mean(explained_peaks_pcts):.1f}%")
            print(f"Median explained peaks: {np.median(explained_peaks_pcts):.1f}%")

        explained_intensity_pcts = [s['explained_intensity_pct'] for s in stats_list if s['explained_intensity_pct'] > 0]
        if explained_intensity_pcts:
            print(f"Mean explained intensity: {np.mean(explained_intensity_pcts):.1f}%")
            print(f"Median explained intensity: {np.median(explained_intensity_pcts):.1f}%")

    print("\nAnalysis complete!")


if __name__ == '__main__':
    main()
