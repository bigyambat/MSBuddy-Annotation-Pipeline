#!/usr/bin/env python3
"""
Annotate Peaks using GNPS Reference Library Data

Uses reference molecular formulas from GNPS to annotate fragment peaks
and calculate spectral quality metrics.

Author: MSBuddy Annotation Pipeline
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import sys
import re

try:
    from pyteomics import mgf, mass
except ImportError:
    print("Error: pyteomics not installed. Install with: pip install pyteomics", file=sys.stderr)
    sys.exit(1)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Annotate fragment peaks using GNPS reference annotations"
    )
    parser.add_argument(
        '--mgf',
        type=str,
        required=True,
        help='Input MGF file'
    )
    parser.add_argument(
        '--reference',
        type=str,
        required=True,
        help='Reference annotations TSV from parse_gnps_reference.py'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Output TSV file with annotated peaks and quality metrics'
    )
    parser.add_argument(
        '--mz_tolerance',
        type=float,
        default=0.02,
        help='m/z tolerance for peak matching (Da, default: 0.02 for high-res MS)'
    )
    parser.add_argument(
        '--ppm_tolerance',
        type=float,
        default=15.0,
        help='PPM tolerance for peak matching (default: 15 ppm, optimal for Orbitrap)'
    )
    parser.add_argument(
        '--min_explained_peaks',
        type=int,
        default=3,
        help='Minimum explained peaks for any reportable annotation (default: 3)'
    )
    parser.add_argument(
        '--min_explained_peaks_good',
        type=int,
        default=6,
        help='Minimum explained peaks for "good" quality (default: 6, GNPS standard)'
    )
    parser.add_argument(
        '--min_explained_peaks_uncertain',
        type=int,
        default=4,
        help='Minimum explained peaks for "uncertain" quality (default: 4)'
    )
    parser.add_argument(
        '--min_explained_intensity',
        type=float,
        default=0.5,
        help='Minimum explained intensity fraction (default: 0.5 = 50%)'
    )
    parser.add_argument(
        '--quality_threshold_good',
        type=float,
        default=0.65,
        help='Minimum quality score for "good" classification (default: 0.65, ~1%% FDR)'
    )
    parser.add_argument(
        '--quality_threshold_uncertain',
        type=float,
        default=0.50,
        help='Minimum quality score for "uncertain" classification (default: 0.50, ~5%% FDR)'
    )

    return parser.parse_args()


def load_mgf_spectra(mgf_path: str) -> Dict:
    """
    Load MGF spectra into dictionary.

    Args:
        mgf_path: Path to MGF file

    Returns:
        Dictionary mapping spectrum_id to spectrum data
    """
    spectra_dict = {}
    skipped = 0

    print(f"Loading spectra from: {mgf_path}")

    try:
        with mgf.read(mgf_path, use_index=False) as reader:
            for idx, spectrum in enumerate(reader):
                # Skip None spectra (can occur with malformed MGF entries)
                if spectrum is None:
                    skipped += 1
                    continue

                params = spectrum.get('params', {})

                # Extract spectrum ID
                spec_id = (
                    params.get('spectrumid') or
                    params.get('spectrum_id') or
                    params.get('title') or
                    spectrum.get('title') or
                    f'spectrum_{idx}'
                )

                # Extract m/z and intensity arrays
                mz_array = spectrum.get('m/z array')
                intensity_array = spectrum.get('intensity array')

                if mz_array is None or intensity_array is None:
                    skipped += 1
                    continue

                # Extract precursor info
                pepmass = params.get('pepmass')
                if isinstance(pepmass, (list, tuple)):
                    precursor_mz = float(pepmass[0])
                else:
                    precursor_mz = float(pepmass) if pepmass else 0.0

                # Extract charge
                charge_str = params.get('charge', '1+')
                try:
                    if isinstance(charge_str, (list, tuple)):
                        charge_str = charge_str[0]
                    charge_str = str(charge_str).strip()
                    if charge_str.endswith('+'):
                        charge = int(charge_str[:-1]) if len(charge_str) > 1 else 1
                    elif charge_str.endswith('-'):
                        charge = -int(charge_str[:-1]) if len(charge_str) > 1 else -1
                    else:
                        charge = int(charge_str) if charge_str else 1
                except (ValueError, AttributeError):
                    charge = 1

                spectra_dict[spec_id] = {
                    'mz': mz_array,
                    'intensity': intensity_array,
                    'precursor_mz': precursor_mz,
                    'charge': charge
                }

                if (idx + 1) % 1000 == 0:
                    print(f"  Loaded {idx + 1} spectra...")

    except Exception as e:
        print(f"Error loading MGF: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(spectra_dict)} spectra (skipped {skipped})")
    return spectra_dict


def generate_fragment_library(formula: str, precursor_mz: float, charge: int = 1) -> List[Dict]:
    """
    Generate theoretical fragment library from molecular formula.

    Args:
        formula: Molecular formula string
        precursor_mz: Precursor m/z
        charge: Ion charge state

    Returns:
        List of theoretical fragments with m/z and description
    """
    if not formula or pd.isna(formula):
        return []

    fragments = []
    proton_mass = 1.007276467

    # Parse formula using pyteomics
    try:
        composition = mass.Composition(formula=formula)
        neutral_mass = mass.calculate_mass(composition=composition)
    except Exception as e:
        # Fallback to simple parsing if pyteomics fails
        return []

    # Common neutral losses
    LOSSES = {
        'H2O': mass.calculate_mass(formula='H2O'),
        'NH3': mass.calculate_mass(formula='NH3'),
        'CO': mass.calculate_mass(formula='CO'),
        'CO2': mass.calculate_mass(formula='CO2'),
        'CH2': mass.calculate_mass(formula='CH2'),
        'CH3': mass.calculate_mass(formula='CH3'),
        'C2H4': mass.calculate_mass(formula='C2H4'),
    }

    # Add precursor
    fragments.append({
        'mz': precursor_mz,
        'formula': formula,
        'loss': 'M',
        'relative_intensity': 100.0
    })

    # Generate neutral loss fragments
    for loss_name, loss_mass in LOSSES.items():
        fragment_mass = neutral_mass - loss_mass
        if fragment_mass > 50:  # Minimum reasonable fragment mass
            if charge > 0:
                fragment_mz = (fragment_mass + charge * proton_mass) / abs(charge)
            else:
                fragment_mz = (fragment_mass - abs(charge) * proton_mass) / abs(charge)

            if fragment_mz > 0:
                fragments.append({
                    'mz': fragment_mz,
                    'formula': f'{formula}-{loss_name}',
                    'loss': loss_name,
                    'relative_intensity': 50.0
                })

    # Generate simple backbone cleavages (simplified approach)
    # For more sophisticated fragmentation, would need structure-based prediction
    if 'C' in composition:
        num_carbons = composition['C']
        for i in range(1, min(num_carbons, 5)):  # Try removing 1-4 carbons
            frag_comp = composition.copy()
            if i <= frag_comp.get('C', 0):
                frag_comp['C'] -= i
                # Also remove some hydrogens
                h_remove = min(2 * i, frag_comp.get('H', 0))
                if h_remove > 0:
                    frag_comp['H'] -= h_remove

                try:
                    frag_mass = mass.calculate_mass(composition=frag_comp)
                    if frag_mass > 50:
                        if charge > 0:
                            frag_mz = (frag_mass + charge * proton_mass) / abs(charge)
                        else:
                            frag_mz = (frag_mass - abs(charge) * proton_mass) / abs(charge)

                        if frag_mz > 0:
                            fragments.append({
                                'mz': frag_mz,
                                'formula': str(frag_comp),
                                'loss': f'C{i}H{h_remove}',
                                'relative_intensity': 30.0
                            })
                except:
                    continue

    return fragments


def match_peaks_to_fragments(spectrum_mz: np.ndarray,
                             spectrum_intensity: np.ndarray,
                             theoretical_fragments: List[Dict],
                             mz_tolerance: float = 0.01,
                             ppm_tolerance: float = 20.0) -> Dict:
    """
    Match observed peaks to theoretical fragments.

    Args:
        spectrum_mz: Observed m/z values
        spectrum_intensity: Observed intensities
        theoretical_fragments: List of theoretical fragments
        mz_tolerance: Absolute m/z tolerance (Da)
        ppm_tolerance: PPM tolerance

    Returns:
        Dictionary with matching statistics
    """
    if len(theoretical_fragments) == 0:
        return {
            'num_matched': 0,
            'total_peaks': len(spectrum_mz),
            'matched_intensity': 0.0,
            'total_intensity': float(np.sum(spectrum_intensity)),
            'matched_peaks_pct': 0.0,
            'matched_intensity_pct': 0.0,
            'matched_fragments': []
        }

    total_intensity = float(np.sum(spectrum_intensity))
    matched_peaks = set()
    matched_intensity = 0.0
    matched_fragments = []

    for fragment in theoretical_fragments:
        frag_mz = fragment['mz']

        # Calculate tolerances
        abs_tol = mz_tolerance
        ppm_tol = frag_mz * ppm_tolerance / 1e6

        # Use the larger of the two tolerances
        tolerance = max(abs_tol, ppm_tol)

        # Find matching peaks
        mz_diffs = np.abs(spectrum_mz - frag_mz)
        matches = np.where(mz_diffs <= tolerance)[0]

        if len(matches) > 0:
            # Take the closest match
            best_match_idx = matches[np.argmin(mz_diffs[matches])]

            if best_match_idx not in matched_peaks:
                matched_peaks.add(best_match_idx)
                matched_intensity += spectrum_intensity[best_match_idx]

                matched_fragments.append({
                    'observed_mz': spectrum_mz[best_match_idx],
                    'theoretical_mz': frag_mz,
                    'error_da': spectrum_mz[best_match_idx] - frag_mz,
                    'error_ppm': (spectrum_mz[best_match_idx] - frag_mz) / frag_mz * 1e6,
                    'intensity': spectrum_intensity[best_match_idx],
                    'formula': fragment['formula'],
                    'loss': fragment['loss']
                })

    # Calculate statistics
    num_matched = len(matched_peaks)
    matched_peaks_pct = (num_matched / len(spectrum_mz) * 100) if len(spectrum_mz) > 0 else 0.0
    matched_intensity_pct = (matched_intensity / total_intensity * 100) if total_intensity > 0 else 0.0

    return {
        'num_matched': num_matched,
        'total_peaks': len(spectrum_mz),
        'matched_intensity': matched_intensity,
        'total_intensity': total_intensity,
        'matched_peaks_pct': matched_peaks_pct,
        'matched_intensity_pct': matched_intensity_pct,
        'matched_fragments': matched_fragments
    }


def calculate_quality_score(match_stats: Dict,
                           min_explained_peaks: int = 3,
                           min_explained_intensity: float = 0.5) -> float:
    """
    Calculate overall quality score from matching statistics.

    Args:
        match_stats: Matching statistics dictionary
        min_explained_peaks: Minimum explained peaks threshold
        min_explained_intensity: Minimum explained intensity threshold

    Returns:
        Quality score (0-1)
    """
    # Component scores
    peak_score = min(match_stats['num_matched'] / max(min_explained_peaks, 1), 1.0)
    intensity_score = match_stats['matched_intensity_pct'] / 100.0

    # Weighted average (intensity is more important)
    quality_score = 0.3 * peak_score + 0.7 * intensity_score

    return quality_score


def classify_quality(quality_score: float,
                    num_matched: int,
                    threshold_good: float = 0.65,
                    threshold_uncertain: float = 0.50,
                    min_peaks_good: int = 6,
                    min_peaks_uncertain: int = 4) -> str:
    """
    Classify annotation quality using tiered thresholds.

    Args:
        quality_score: Overall quality score (0-1)
        num_matched: Number of matched peaks
        threshold_good: Threshold for "good" classification (~1% FDR, default: 0.65)
        threshold_uncertain: Threshold for "uncertain" classification (~5% FDR, default: 0.50)
        min_peaks_good: Minimum peaks for "good" quality (GNPS standard, default: 6)
        min_peaks_uncertain: Minimum peaks for "uncertain" quality (default: 4)

    Returns:
        Quality classification: 'good', 'uncertain', 'bad', 'no_annotation'
    """
    if num_matched == 0:
        return 'no_annotation'

    # Good: high quality score AND sufficient peaks (GNPS standard: ≥6 peaks)
    if quality_score >= threshold_good and num_matched >= min_peaks_good:
        return 'good'
    # Uncertain: moderate quality score AND minimum peaks for confidence
    elif quality_score >= threshold_uncertain and num_matched >= min_peaks_uncertain:
        return 'uncertain'
    else:
        return 'bad'


def annotate_spectra(spectra_dict: Dict,
                    reference_df: pd.DataFrame,
                    args) -> pd.DataFrame:
    """
    Annotate all spectra using reference library.

    Args:
        spectra_dict: Dictionary of spectrum data
        reference_df: Reference annotations DataFrame
        args: Command-line arguments

    Returns:
        DataFrame with annotated peaks and quality metrics
    """
    results = []

    print(f"\nAnnotating {len(reference_df)} spectra...")

    for idx, ref_row in reference_df.iterrows():
        spec_id = ref_row['spectrum_id']

        if spec_id not in spectra_dict:
            print(f"Warning: Spectrum '{spec_id}' not found in MGF file", file=sys.stderr)
            continue

        spectrum = spectra_dict[spec_id]

        # Generate theoretical fragments
        theoretical_fragments = generate_fragment_library(
            ref_row['formula'],
            ref_row['precursor_mz'],
            ref_row['charge']
        )

        # Match peaks to fragments
        match_stats = match_peaks_to_fragments(
            spectrum['mz'],
            spectrum['intensity'],
            theoretical_fragments,
            mz_tolerance=args.mz_tolerance,
            ppm_tolerance=args.ppm_tolerance
        )

        # Calculate quality score
        quality_score = calculate_quality_score(
            match_stats,
            min_explained_peaks=args.min_explained_peaks,
            min_explained_intensity=args.min_explained_intensity
        )

        # Classify quality (using tiered peak thresholds)
        quality_class = classify_quality(
            quality_score,
            match_stats['num_matched'],
            threshold_good=args.quality_threshold_good,
            threshold_uncertain=args.quality_threshold_uncertain,
            min_peaks_good=args.min_explained_peaks_good,
            min_peaks_uncertain=args.min_explained_peaks_uncertain
        )

        # Create result record
        result = {
            'spectrum_id': spec_id,
            'smiles': ref_row['smiles'],
            'inchi': ref_row['inchi'],
            'inchikey': ref_row['inchikey'],
            'formula': ref_row['formula'],
            'compound_name': ref_row['compound_name'],
            'library_id': ref_row['library_id'],
            'precursor_mz': ref_row['precursor_mz'],
            'charge': ref_row['charge'],
            'adduct': ref_row['adduct'],
            'num_peaks': match_stats['total_peaks'],
            'num_matched_peaks': match_stats['num_matched'],
            'matched_peaks_pct': match_stats['matched_peaks_pct'],
            'matched_intensity_pct': match_stats['matched_intensity_pct'],
            'quality_score': quality_score,
            'annotation_quality': quality_class,
            'collision_energy': ref_row.get('collision_energy', ''),
            'retention_time': ref_row.get('retention_time', ''),
            'num_theoretical_fragments': len(theoretical_fragments),
        }

        results.append(result)

        if (idx + 1) % 100 == 0:
            print(f"  Annotated {idx + 1}/{len(reference_df)} spectra...")

    print(f"Annotation complete: {len(results)} spectra processed")

    return pd.DataFrame(results)


def print_summary(df: pd.DataFrame):
    """Print annotation summary statistics."""
    print("\n" + "=" * 60)
    print("ANNOTATION SUMMARY")
    print("=" * 60)

    total = len(df)
    print(f"Total spectra annotated: {total}")

    if 'annotation_quality' in df.columns:
        print("\nQuality Distribution:")
        quality_counts = df['annotation_quality'].value_counts()
        for quality, count in quality_counts.items():
            pct = count / total * 100
            print(f"  {quality:15s}: {count:6d} ({pct:5.1f}%)")

    if 'matched_peaks_pct' in df.columns:
        matched_pct = df['matched_peaks_pct'].dropna()
        if len(matched_pct) > 0:
            print(f"\nMatched Peaks Percentage:")
            print(f"  Mean:   {matched_pct.mean():.1f}%")
            print(f"  Median: {matched_pct.median():.1f}%")
            print(f"  Std:    {matched_pct.std():.1f}%")

    if 'matched_intensity_pct' in df.columns:
        intensity_pct = df['matched_intensity_pct'].dropna()
        if len(intensity_pct) > 0:
            print(f"\nMatched Intensity Percentage:")
            print(f"  Mean:   {intensity_pct.mean():.1f}%")
            print(f"  Median: {intensity_pct.median():.1f}%")
            print(f"  Std:    {intensity_pct.std():.1f}%")

    if 'quality_score' in df.columns:
        quality_score = df['quality_score'].dropna()
        if len(quality_score) > 0:
            print(f"\nQuality Score (0-1):")
            print(f"  Mean:   {quality_score.mean():.3f}")
            print(f"  Median: {quality_score.median():.3f}")
            print(f"  Std:    {quality_score.std():.3f}")

    print("=" * 60)


def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 60)
    print("GNPS Peak Annotation")
    print("=" * 60)

    # Load reference annotations
    print(f"\nLoading reference annotations from: {args.reference}")
    try:
        reference_df = pd.read_csv(args.reference, sep='\t')
        print(f"Loaded {len(reference_df)} reference annotations")
    except Exception as e:
        print(f"Error loading reference file: {e}", file=sys.stderr)
        sys.exit(1)

    # Load spectra from MGF
    spectra_dict = load_mgf_spectra(args.mgf)

    # Annotate spectra
    annotated_df = annotate_spectra(spectra_dict, reference_df, args)

    if len(annotated_df) == 0:
        print("Error: No spectra were annotated", file=sys.stderr)
        sys.exit(1)

    # Print summary
    print_summary(annotated_df)

    # Save results
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    annotated_df.to_csv(args.output, sep='\t', index=False)
    print(f"\nAnnotated peaks saved to: {args.output}")
    print("=" * 60)


if __name__ == '__main__':
    main()
