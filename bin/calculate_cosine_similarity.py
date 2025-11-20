#!/usr/bin/env python3
"""
Calculate Cosine Similarity Between Mass Spectra

Computes pairwise spectral similarity using modified cosine score,
identifies spectral clusters, and finds most similar neighbors.

Author: MSBuddy Annotation Pipeline
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import sys
import os
from collections import defaultdict
from multiprocessing import Pool, cpu_count

try:
    from pyteomics import mgf
except ImportError:
    print("Error: pyteomics not installed. Install with: pip install pyteomics", file=sys.stderr)
    sys.exit(1)

# Optional: LSH for fast approximate nearest neighbor filtering
try:
    from datasketch import MinHash, MinHashLSH
    DATASKETCH_AVAILABLE = True
except ImportError:
    DATASKETCH_AVAILABLE = False
    print("Warning: datasketch not installed. LSH optimization disabled.", file=sys.stderr)

# Optional: joblib for parallel processing
try:
    from joblib import Parallel, delayed
    JOBLIB_AVAILABLE = True
except ImportError:
    JOBLIB_AVAILABLE = False
    print("Warning: joblib not installed. Using multiprocessing instead.", file=sys.stderr)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate cosine similarity between MS/MS spectra"
    )
    parser.add_argument(
        '--mgf',
        type=str,
        required=True,
        help='Input MGF file'
    )
    parser.add_argument(
        '--annotations',
        type=str,
        required=True,
        help='Annotated peaks TSV from annotate_peaks_gnps.py'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Output TSV file with similarity scores'
    )
    parser.add_argument(
        '--output_matrix',
        type=str,
        default=None,
        help='Optional: Output full similarity matrix as CSV'
    )
    parser.add_argument(
        '--mz_tolerance',
        type=float,
        default=0.01,
        help='m/z tolerance for peak matching (Da, default: 0.01)'
    )
    parser.add_argument(
        '--min_peaks',
        type=int,
        default=3,
        help='Minimum matching peaks for similarity calculation (default: 3)'
    )
    parser.add_argument(
        '--top_n',
        type=int,
        default=10,
        help='Number of top similar spectra to report per spectrum (default: 10)'
    )
    parser.add_argument(
        '--min_similarity',
        type=float,
        default=0.5,
        help='Minimum similarity score to report (default: 0.5)'
    )
    parser.add_argument(
        '--intensity_power',
        type=float,
        default=0.5,
        help='Power to scale intensities (default: 0.5, use 1.0 for no scaling)'
    )
    parser.add_argument(
        '--use_lsh',
        action='store_true',
        default=False,
        help='Use LSH for fast approximate filtering (recommended for >1000 spectra)'
    )
    parser.add_argument(
        '--lsh_threshold',
        type=float,
        default=0.4,
        help='LSH similarity threshold for candidate filtering (default: 0.4)'
    )
    parser.add_argument(
        '--lsh_num_perm',
        type=int,
        default=128,
        help='Number of permutations for MinHash (default: 128)'
    )
    parser.add_argument(
        '--num_workers',
        type=int,
        default=0,
        help='Number of parallel workers (default: 0 = auto-detect CPU count)'
    )
    parser.add_argument(
        '--chunk_size',
        type=int,
        default=1000,
        help='Chunk size for parallel processing (default: 1000)'
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

                spectra_dict[spec_id] = {
                    'mz': mz_array,
                    'intensity': intensity_array,
                    'precursor_mz': precursor_mz
                }

                if (idx + 1) % 1000 == 0:
                    print(f"  Loaded {idx + 1} spectra...")

    except Exception as e:
        print(f"Error loading MGF: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(spectra_dict)} spectra (skipped {skipped})")
    return spectra_dict


def normalize_spectrum(mz: np.ndarray, intensity: np.ndarray, power: float = 0.5) -> Tuple[np.ndarray, np.ndarray]:
    """
    Normalize spectrum intensities.

    Args:
        mz: m/z array
        intensity: Intensity array
        power: Power to scale intensities (0.5 = sqrt scaling)

    Returns:
        Tuple of (mz, normalized_intensity)
    """
    if len(intensity) == 0:
        return mz, intensity

    # Apply power scaling
    scaled_intensity = np.power(intensity, power)

    # Normalize to unit vector
    norm = np.linalg.norm(scaled_intensity)
    if norm > 0:
        normalized_intensity = scaled_intensity / norm
    else:
        normalized_intensity = scaled_intensity

    return mz, normalized_intensity


def spectrum_to_minhash(mz: np.ndarray, intensity: np.ndarray,
                        num_perm: int = 128, mz_bin_size: float = 0.01):
    """
    Convert a spectrum to a MinHash signature for LSH.

    Args:
        mz: m/z array
        intensity: Intensity array
        num_perm: Number of permutations for MinHash
        mz_bin_size: Bin size for discretizing m/z values

    Returns:
        MinHash signature
    """
    m = MinHash(num_perm=num_perm)

    if len(mz) == 0:
        return m

    # Weight peaks by intensity (top peaks get more hash entries)
    max_intensity = np.max(intensity) if len(intensity) > 0 else 1.0

    for i, (mz_val, int_val) in enumerate(zip(mz, intensity)):
        # Discretize m/z to bin
        mz_bin = int(mz_val / mz_bin_size)

        # Add hash for this peak
        m.update(str(mz_bin).encode('utf8'))

        # Add extra hashes for high-intensity peaks (improves accuracy)
        relative_int = int_val / max_intensity if max_intensity > 0 else 0
        if relative_int > 0.5:
            m.update(f"{mz_bin}_high".encode('utf8'))
        if relative_int > 0.8:
            m.update(f"{mz_bin}_top".encode('utf8'))

    return m


def build_lsh_index(normalized_spectra: Dict, num_perm: int = 128,
                    threshold: float = 0.4) -> Tuple:
    """
    Build LSH index for fast candidate pair filtering.

    Args:
        normalized_spectra: Dictionary of normalized spectrum data
        num_perm: Number of permutations for MinHash
        threshold: Jaccard similarity threshold for LSH

    Returns:
        Tuple of (LSH index, spectrum_id -> MinHash mapping)
    """
    if not DATASKETCH_AVAILABLE:
        return None, {}

    print(f"  Building LSH index (threshold={threshold}, num_perm={num_perm})...")

    lsh = MinHashLSH(threshold=threshold, num_perm=num_perm)
    minhashes = {}

    for spec_id, spec in normalized_spectra.items():
        mh = spectrum_to_minhash(spec['mz'], spec['intensity'], num_perm=num_perm)
        minhashes[spec_id] = mh
        lsh.insert(spec_id, mh)

    return lsh, minhashes


def get_lsh_candidate_pairs(lsh, minhashes: Dict,
                            spectrum_ids: List[str]) -> List[Tuple[str, str]]:
    """
    Get candidate pairs from LSH index.

    Args:
        lsh: MinHashLSH index
        minhashes: Dictionary of MinHash signatures
        spectrum_ids: List of spectrum IDs

    Returns:
        List of candidate pairs (spec_id1, spec_id2)
    """
    print("  Querying LSH for candidate pairs...")

    candidates = set()

    for spec_id in spectrum_ids:
        # Query LSH for similar spectra
        similar = lsh.query(minhashes[spec_id])

        for match_id in similar:
            if match_id != spec_id:
                # Store as sorted tuple to avoid duplicates
                pair = tuple(sorted([spec_id, match_id]))
                candidates.add(pair)

    candidate_list = list(candidates)
    total_possible = len(spectrum_ids) * (len(spectrum_ids) - 1) // 2
    reduction = (1 - len(candidate_list) / total_possible) * 100 if total_possible > 0 else 0

    print(f"  LSH filtered to {len(candidate_list)} candidate pairs ({reduction:.1f}% reduction)")

    return candidate_list


def compute_similarity_for_pair(args_tuple):
    """
    Compute cosine similarity for a single pair (for parallel processing).

    Args:
        args_tuple: Tuple of (spec_id1, spec_id2, spec1, spec2, mz_tolerance, min_peaks, min_similarity)

    Returns:
        Dictionary with similarity result or None
    """
    spec_id1, spec_id2, spec1, spec2, mz_tolerance, min_peaks, min_similarity = args_tuple

    cosine_score, num_matches = calculate_cosine_similarity(
        spec1['mz'], spec1['intensity'],
        spec2['mz'], spec2['intensity'],
        mz_tolerance=mz_tolerance,
        min_peaks=min_peaks
    )

    if cosine_score >= min_similarity:
        return {
            'spectrum_id_1': spec_id1,
            'spectrum_id_2': spec_id2,
            'cosine_score': cosine_score,
            'num_matching_peaks': num_matches,
            'precursor_mz_1': spec1['precursor_mz'],
            'precursor_mz_2': spec2['precursor_mz'],
            'precursor_mz_diff': abs(spec1['precursor_mz'] - spec2['precursor_mz'])
        }
    return None


def calculate_cosine_similarity(mz1: np.ndarray,
                                intensity1: np.ndarray,
                                mz2: np.ndarray,
                                intensity2: np.ndarray,
                                mz_tolerance: float = 0.01,
                                min_peaks: int = 3) -> Tuple[float, int]:
    """
    Calculate modified cosine similarity between two spectra.

    Args:
        mz1: m/z array of spectrum 1
        intensity1: Intensity array of spectrum 1 (normalized)
        mz2: m/z array of spectrum 2
        intensity2: Intensity array of spectrum 2 (normalized)
        mz_tolerance: m/z tolerance for peak matching (Da)
        min_peaks: Minimum matching peaks required

    Returns:
        Tuple of (cosine_score, num_matching_peaks)
    """
    if len(mz1) == 0 or len(mz2) == 0:
        return 0.0, 0

    # Find matching peaks
    matched_intensity1 = []
    matched_intensity2 = []

    # Use a simple approach: for each peak in spectrum 1, find closest in spectrum 2
    for i, m1 in enumerate(mz1):
        diffs = np.abs(mz2 - m1)
        min_diff_idx = np.argmin(diffs)

        if diffs[min_diff_idx] <= mz_tolerance:
            matched_intensity1.append(intensity1[i])
            matched_intensity2.append(intensity2[min_diff_idx])

    num_matches = len(matched_intensity1)

    if num_matches < min_peaks:
        return 0.0, num_matches

    # Calculate cosine similarity
    matched_intensity1 = np.array(matched_intensity1)
    matched_intensity2 = np.array(matched_intensity2)

    # Cosine similarity = dot product (since vectors are normalized)
    cosine_score = np.dot(matched_intensity1, matched_intensity2)

    # Clip to [0, 1] range (numerical errors can cause slight exceeding)
    cosine_score = np.clip(cosine_score, 0.0, 1.0)

    return float(cosine_score), num_matches


def calculate_pairwise_similarities(spectra_dict: Dict,
                                   args) -> pd.DataFrame:
    """
    Calculate pairwise cosine similarities between all spectra.

    Supports LSH pre-filtering and parallel processing for large datasets.

    Args:
        spectra_dict: Dictionary of spectrum data
        args: Command-line arguments

    Returns:
        DataFrame with pairwise similarities
    """
    spectrum_ids = list(spectra_dict.keys())
    n_spectra = len(spectrum_ids)
    total_comparisons = n_spectra * (n_spectra - 1) // 2

    print(f"\nCalculating pairwise similarities for {n_spectra} spectra...")
    print(f"Total possible comparisons: {total_comparisons}")

    # Normalize all spectra first
    print("  Normalizing spectra...")
    normalized_spectra = {}
    for spec_id, spec in spectra_dict.items():
        mz_norm, intensity_norm = normalize_spectrum(
            spec['mz'],
            spec['intensity'],
            power=args.intensity_power
        )
        normalized_spectra[spec_id] = {
            'mz': mz_norm,
            'intensity': intensity_norm,
            'precursor_mz': spec['precursor_mz']
        }

    # Determine number of workers
    num_workers = args.num_workers if args.num_workers > 0 else cpu_count()

    # Use LSH for large datasets if enabled
    use_lsh = args.use_lsh and DATASKETCH_AVAILABLE and n_spectra > 100

    if use_lsh:
        print(f"\nUsing LSH optimization (threshold={args.lsh_threshold})...")
        lsh, minhashes = build_lsh_index(
            normalized_spectra,
            num_perm=args.lsh_num_perm,
            threshold=args.lsh_threshold
        )
        candidate_pairs = get_lsh_candidate_pairs(lsh, minhashes, spectrum_ids)
    else:
        # Generate all pairs
        print("  Generating all pairs (no LSH filtering)...")
        candidate_pairs = []
        for i in range(n_spectra):
            for j in range(i + 1, n_spectra):
                candidate_pairs.append((spectrum_ids[i], spectrum_ids[j]))

    num_candidates = len(candidate_pairs)
    print(f"  Processing {num_candidates} candidate pairs...")

    # Decide between parallel and sequential processing
    use_parallel = num_candidates > 1000 and num_workers > 1

    if use_parallel:
        print(f"  Using parallel processing with {num_workers} workers...")

        # Prepare arguments for parallel processing
        pair_args = [
            (
                spec_id1, spec_id2,
                normalized_spectra[spec_id1], normalized_spectra[spec_id2],
                args.mz_tolerance, args.min_peaks, args.min_similarity
            )
            for spec_id1, spec_id2 in candidate_pairs
        ]

        # Use joblib if available, otherwise multiprocessing
        if JOBLIB_AVAILABLE:
            results = Parallel(n_jobs=num_workers, verbose=1)(
                delayed(compute_similarity_for_pair)(pair_arg)
                for pair_arg in pair_args
            )
        else:
            # Use multiprocessing Pool
            with Pool(num_workers) as pool:
                chunk_size = max(1, len(pair_args) // (num_workers * 10))
                results = []
                for i, result in enumerate(pool.imap_unordered(
                    compute_similarity_for_pair, pair_args, chunksize=chunk_size
                )):
                    results.append(result)
                    if (i + 1) % 10000 == 0:
                        progress = (i + 1) / num_candidates * 100
                        print(f"    Progress: {i + 1}/{num_candidates} ({progress:.1f}%)")

        # Filter None results
        similarities = [r for r in results if r is not None]

    else:
        # Sequential processing
        similarities = []
        for i, (spec_id1, spec_id2) in enumerate(candidate_pairs):
            spec1 = normalized_spectra[spec_id1]
            spec2 = normalized_spectra[spec_id2]

            cosine_score, num_matches = calculate_cosine_similarity(
                spec1['mz'],
                spec1['intensity'],
                spec2['mz'],
                spec2['intensity'],
                mz_tolerance=args.mz_tolerance,
                min_peaks=args.min_peaks
            )

            if cosine_score >= args.min_similarity:
                similarities.append({
                    'spectrum_id_1': spec_id1,
                    'spectrum_id_2': spec_id2,
                    'cosine_score': cosine_score,
                    'num_matching_peaks': num_matches,
                    'precursor_mz_1': spec1['precursor_mz'],
                    'precursor_mz_2': spec2['precursor_mz'],
                    'precursor_mz_diff': abs(spec1['precursor_mz'] - spec2['precursor_mz'])
                })

            if (i + 1) % 10000 == 0:
                progress = (i + 1) / num_candidates * 100
                print(f"  Progress: {i + 1}/{num_candidates} ({progress:.1f}%)")

    print(f"Found {len(similarities)} similar pairs (score >= {args.min_similarity})")

    return pd.DataFrame(similarities)


def find_top_neighbors(similarity_df: pd.DataFrame, top_n: int = 10) -> pd.DataFrame:
    """
    Find top N most similar neighbors for each spectrum.

    Args:
        similarity_df: DataFrame with pairwise similarities
        top_n: Number of top neighbors to keep

    Returns:
        DataFrame with top neighbors per spectrum
    """
    if len(similarity_df) == 0:
        return pd.DataFrame()

    print(f"\nFinding top {top_n} neighbors for each spectrum...")

    # Create bidirectional similarity records
    records = []

    for _, row in similarity_df.iterrows():
        # Add both directions
        records.append({
            'query_spectrum': row['spectrum_id_1'],
            'match_spectrum': row['spectrum_id_2'],
            'cosine_score': row['cosine_score'],
            'num_matching_peaks': row['num_matching_peaks'],
            'precursor_mz_diff': row['precursor_mz_diff']
        })
        records.append({
            'query_spectrum': row['spectrum_id_2'],
            'match_spectrum': row['spectrum_id_1'],
            'cosine_score': row['cosine_score'],
            'num_matching_peaks': row['num_matching_peaks'],
            'precursor_mz_diff': row['precursor_mz_diff']
        })

    all_similarities = pd.DataFrame(records)

    # Group by query and take top N
    top_neighbors = []

    for query_spec in all_similarities['query_spectrum'].unique():
        matches = all_similarities[all_similarities['query_spectrum'] == query_spec]
        matches_sorted = matches.sort_values('cosine_score', ascending=False).head(top_n)

        for rank, (_, match_row) in enumerate(matches_sorted.iterrows(), 1):
            top_neighbors.append({
                'query_spectrum': query_spec,
                'match_spectrum': match_row['match_spectrum'],
                'rank': rank,
                'cosine_score': match_row['cosine_score'],
                'num_matching_peaks': match_row['num_matching_peaks'],
                'precursor_mz_diff': match_row['precursor_mz_diff']
            })

    return pd.DataFrame(top_neighbors)


def create_similarity_matrix(similarity_df: pd.DataFrame,
                             spectrum_ids: List[str]) -> np.ndarray:
    """
    Create full similarity matrix from pairwise similarities.

    Args:
        similarity_df: DataFrame with pairwise similarities
        spectrum_ids: List of all spectrum IDs

    Returns:
        NxN similarity matrix
    """
    n = len(spectrum_ids)
    id_to_idx = {spec_id: i for i, spec_id in enumerate(spectrum_ids)}

    # Initialize matrix with zeros
    matrix = np.zeros((n, n))

    # Diagonal is 1.0 (self-similarity)
    np.fill_diagonal(matrix, 1.0)

    # Fill in pairwise similarities (symmetric)
    for _, row in similarity_df.iterrows():
        i = id_to_idx.get(row['spectrum_id_1'])
        j = id_to_idx.get(row['spectrum_id_2'])

        if i is not None and j is not None:
            matrix[i, j] = row['cosine_score']
            matrix[j, i] = row['cosine_score']  # Symmetric

    return matrix


def enhance_annotations_with_similarity(annotations_df: pd.DataFrame,
                                        top_neighbors_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add similarity information to annotations DataFrame.

    Args:
        annotations_df: Annotated peaks DataFrame
        top_neighbors_df: Top neighbors DataFrame

    Returns:
        Enhanced annotations DataFrame
    """
    if len(top_neighbors_df) == 0:
        annotations_df['num_similar_neighbors'] = 0
        annotations_df['max_similarity_score'] = 0.0
        annotations_df['avg_similarity_score'] = 0.0
        return annotations_df

    # Calculate per-spectrum similarity statistics
    similarity_stats = []

    for spec_id in annotations_df['spectrum_id']:
        neighbors = top_neighbors_df[top_neighbors_df['query_spectrum'] == spec_id]

        if len(neighbors) > 0:
            similarity_stats.append({
                'spectrum_id': spec_id,
                'num_similar_neighbors': len(neighbors),
                'max_similarity_score': neighbors['cosine_score'].max(),
                'avg_similarity_score': neighbors['cosine_score'].mean(),
                'top_neighbor': neighbors.iloc[0]['match_spectrum'] if len(neighbors) > 0 else '',
                'top_neighbor_score': neighbors.iloc[0]['cosine_score'] if len(neighbors) > 0 else 0.0
            })
        else:
            similarity_stats.append({
                'spectrum_id': spec_id,
                'num_similar_neighbors': 0,
                'max_similarity_score': 0.0,
                'avg_similarity_score': 0.0,
                'top_neighbor': '',
                'top_neighbor_score': 0.0
            })

    similarity_stats_df = pd.DataFrame(similarity_stats)

    # Merge with annotations
    enhanced_df = annotations_df.merge(
        similarity_stats_df,
        on='spectrum_id',
        how='left'
    )

    return enhanced_df


def print_summary(similarity_df: pd.DataFrame, top_neighbors_df: pd.DataFrame):
    """Print similarity calculation summary."""
    print("\n" + "=" * 60)
    print("SIMILARITY CALCULATION SUMMARY")
    print("=" * 60)

    if len(similarity_df) > 0:
        print(f"Total similar pairs: {len(similarity_df)}")
        print(f"\nCosine Score Statistics:")
        print(f"  Mean:   {similarity_df['cosine_score'].mean():.3f}")
        print(f"  Median: {similarity_df['cosine_score'].median():.3f}")
        print(f"  Min:    {similarity_df['cosine_score'].min():.3f}")
        print(f"  Max:    {similarity_df['cosine_score'].max():.3f}")

        print(f"\nMatching Peaks Statistics:")
        print(f"  Mean:   {similarity_df['num_matching_peaks'].mean():.1f}")
        print(f"  Median: {similarity_df['num_matching_peaks'].median():.0f}")

    if len(top_neighbors_df) > 0:
        unique_queries = top_neighbors_df['query_spectrum'].nunique()
        print(f"\nSpectra with neighbors: {unique_queries}")

        neighbors_per_spec = top_neighbors_df.groupby('query_spectrum').size()
        print(f"Neighbors per spectrum:")
        print(f"  Mean:   {neighbors_per_spec.mean():.1f}")
        print(f"  Median: {neighbors_per_spec.median():.0f}")

    print("=" * 60)


def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 60)
    print("Spectral Cosine Similarity Calculation")
    print("=" * 60)

    # Load annotations
    print(f"\nLoading annotations from: {args.annotations}")
    try:
        annotations_df = pd.read_csv(args.annotations, sep='\t')
        print(f"Loaded {len(annotations_df)} annotations")
    except Exception as e:
        print(f"Error loading annotations: {e}", file=sys.stderr)
        sys.exit(1)

    # Load spectra from MGF
    spectra_dict = load_mgf_spectra(args.mgf)

    # Calculate pairwise similarities
    similarity_df = calculate_pairwise_similarities(spectra_dict, args)

    # Find top neighbors
    top_neighbors_df = find_top_neighbors(similarity_df, top_n=args.top_n)

    # Print summary
    print_summary(similarity_df, top_neighbors_df)

    # Enhance annotations with similarity info
    enhanced_annotations = enhance_annotations_with_similarity(
        annotations_df,
        top_neighbors_df
    )

    # Save results
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save enhanced annotations
    enhanced_annotations.to_csv(args.output, sep='\t', index=False)
    print(f"\nEnhanced annotations saved to: {args.output}")

    # Save top neighbors
    neighbors_output = output_path.parent / f"{output_path.stem}_neighbors.tsv"
    top_neighbors_df.to_csv(neighbors_output, sep='\t', index=False)
    print(f"Top neighbors saved to: {neighbors_output}")

    # Save full similarity matrix if requested
    if args.output_matrix:
        spectrum_ids = list(spectra_dict.keys())
        similarity_matrix = create_similarity_matrix(similarity_df, spectrum_ids)

        matrix_df = pd.DataFrame(
            similarity_matrix,
            index=spectrum_ids,
            columns=spectrum_ids
        )
        matrix_df.to_csv(args.output_matrix, sep='\t')
        print(f"Similarity matrix saved to: {args.output_matrix}")

    print("=" * 60)


if __name__ == '__main__':
    main()
