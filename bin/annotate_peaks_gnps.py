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

# RDKit for structure-aware fragmentation
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Fragments
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available. Structure-aware fragmentation disabled.", file=sys.stderr)


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


def generate_structure_fragments(smiles: str, precursor_mz: float, charge: int = 1) -> List[Dict]:
    """
    Generate structure-aware fragments using RDKit functional group analysis.

    Identifies functional groups and generates fragments based on common
    cleavage patterns for each functional group type.

    Args:
        smiles: SMILES string of the molecule
        precursor_mz: Precursor m/z
        charge: Ion charge state

    Returns:
        List of theoretical fragments with m/z and description
    """
    if not RDKIT_AVAILABLE or not smiles:
        return []

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []
    except:
        return []

    fragments = []
    proton_mass = 1.007276467
    neutral_mass = Descriptors.ExactMolWt(mol)

    # Define functional group SMARTS patterns and their characteristic losses
    FUNCTIONAL_GROUP_CLEAVAGES = {
        # Esters: lose alkoxy group or acyl group
        'ester': {
            'smarts': '[CX3](=O)[OX2][#6]',
            'losses': [
                ('HCOO', 44.997655),      # Formate loss
                ('CH3COO', 59.013305),    # Acetate loss
                ('CO2', 43.989829),       # CO2 from decarboxylation
                ('OCH3', 31.018175),      # Methoxy loss
                ('OC2H5', 45.033825),     # Ethoxy loss
            ]
        },
        # Amides: lose NH3, CONH2, or alkyl amide
        'amide': {
            'smarts': '[CX3](=O)[NX3]',
            'losses': [
                ('NH3', 17.026549),       # Ammonia loss
                ('CONH2', 44.013663),     # Carboxamide loss
                ('HCONH2', 45.021464),    # Formamide loss
                ('CH3CONH2', 59.037114),  # Acetamide loss
            ]
        },
        # Carboxylic acids: lose H2O, CO2, COOH
        'carboxylic_acid': {
            'smarts': '[CX3](=O)[OX2H1]',
            'losses': [
                ('H2O', 18.010565),       # Water loss
                ('CO2', 43.989829),       # Decarboxylation
                ('COOH', 44.997655),      # Carboxyl radical loss
                ('CO', 27.994915),        # CO loss
            ]
        },
        # Alcohols: lose H2O, or alkyl alcohol
        'alcohol': {
            'smarts': '[OX2H][CX4]',
            'losses': [
                ('H2O', 18.010565),       # Water loss
                ('CH3OH', 32.026215),     # Methanol loss
                ('C2H5OH', 46.041865),    # Ethanol loss
            ]
        },
        # Phenols: lose CO, CHO
        'phenol': {
            'smarts': 'c[OX2H]',
            'losses': [
                ('CO', 27.994915),        # CO loss from ring
                ('CHO', 29.002740),       # Formyl loss
                ('H2O', 18.010565),       # Water loss
            ]
        },
        # Ethers: C-O-C cleavage
        'ether': {
            'smarts': '[#6][OX2][#6]',
            'losses': [
                ('CH2O', 30.010565),      # Formaldehyde loss
                ('C2H4O', 44.026215),     # Acetaldehyde loss
                ('OCH3', 31.018175),      # Methoxy radical
            ]
        },
        # Ketones: lose CO, alkyl groups
        'ketone': {
            'smarts': '[#6][CX3](=O)[#6]',
            'losses': [
                ('CO', 27.994915),        # CO loss
                ('CH3CO', 43.018390),     # Acetyl loss
                ('C2H5CO', 57.034040),    # Propionyl loss
            ]
        },
        # Aldehydes: lose CHO, CO
        'aldehyde': {
            'smarts': '[CX3H1](=O)[#6]',
            'losses': [
                ('CHO', 29.002740),       # Formyl loss
                ('CO', 27.994915),        # CO loss
                ('H2CO', 30.010565),      # Formaldehyde
            ]
        },
        # Primary amines: lose NH3, CH2NH2
        'primary_amine': {
            'smarts': '[NX3H2][CX4]',
            'losses': [
                ('NH3', 17.026549),       # Ammonia loss
                ('CH2NH2', 30.034374),    # Methylamine radical
                ('CH5N', 31.042199),      # Methylamine
            ]
        },
        # Secondary amines
        'secondary_amine': {
            'smarts': '[NX3H1]([#6])[#6]',
            'losses': [
                ('NH3', 17.026549),       # Ammonia loss (with H transfer)
                ('CH3NH', 29.026549),     # Methylenimine
                ('C2H5NH', 43.042199),    # Ethylenimine
            ]
        },
        # Thiols: lose H2S, SH
        'thiol': {
            'smarts': '[SX2H][#6]',
            'losses': [
                ('H2S', 33.987721),       # Hydrogen sulfide
                ('SH', 32.979846),        # Thiol radical
                ('CH3SH', 48.003371),     # Methanethiol
            ]
        },
        # Sulfides/Thioethers
        'sulfide': {
            'smarts': '[#6][SX2][#6]',
            'losses': [
                ('CH2S', 45.987721),      # Thioformaldehyde
                ('SCH3', 46.995546),      # Methylthio radical
            ]
        },
        # Nitriles: lose HCN
        'nitrile': {
            'smarts': '[CX2]#N',
            'losses': [
                ('HCN', 27.010899),       # Hydrogen cyanide
                ('CN', 26.003074),        # Cyano radical
            ]
        },
        # Nitro groups: lose NO2, NO
        'nitro': {
            'smarts': '[NX3](=O)=O',
            'losses': [
                ('NO2', 45.992904),       # Nitrogen dioxide
                ('NO', 29.997989),        # Nitric oxide
                ('HNO2', 46.005479),      # Nitrous acid
            ]
        },
        # Phosphate esters
        'phosphate': {
            'smarts': '[PX4](=O)([OX2])([OX2])[OX2]',
            'losses': [
                ('H3PO4', 97.976896),     # Phosphoric acid
                ('HPO3', 79.966331),      # Metaphosphoric acid
                ('H2PO4', 96.969021),     # Dihydrogen phosphate
            ]
        },
        # Sulfates
        'sulfate': {
            'smarts': '[SX4](=O)(=O)([OX2])[OX2]',
            'losses': [
                ('SO3', 79.956815),       # Sulfur trioxide
                ('H2SO4', 97.967379),     # Sulfuric acid
                ('HSO3', 80.964640),      # Bisulfite
            ]
        },
        # Glycosidic bonds (sugars)
        'glycoside': {
            'smarts': '[OX2]([#6])[CX4]1[OX2][CX4][CX4][CX4][CX4]1',
            'losses': [
                ('C6H10O5', 162.052824),  # Hexose (glucose/galactose)
                ('C5H8O4', 132.042259),   # Pentose (xylose/arabinose)
                ('C6H10O4', 146.057909),  # Deoxyhexose (fucose/rhamnose)
            ]
        },
        # Aromatic rings - characteristic fragments
        'benzene': {
            'smarts': 'c1ccccc1',
            'losses': [
                ('C2H2', 26.015650),      # Acetylene (retro-Diels-Alder)
                ('C4H2', 50.015650),      # Diacetylene
                ('CO', 27.994915),        # CO from substituted rings
            ]
        },
        # Halides
        'chloride': {
            'smarts': '[Cl][#6]',
            'losses': [
                ('HCl', 35.976678),       # Hydrogen chloride
                ('Cl', 34.968853),        # Chlorine radical
            ]
        },
        'bromide': {
            'smarts': '[Br][#6]',
            'losses': [
                ('HBr', 79.926160),       # Hydrogen bromide
                ('Br', 78.918336),        # Bromine radical
            ]
        },
        'fluoride': {
            'smarts': '[F][#6]',
            'losses': [
                ('HF', 20.006229),        # Hydrogen fluoride
                ('F', 18.998403),         # Fluorine radical
            ]
        },
    }

    # Find functional groups present in the molecule
    found_groups = []
    for group_name, group_info in FUNCTIONAL_GROUP_CLEAVAGES.items():
        pattern = Chem.MolFromSmarts(group_info['smarts'])
        if pattern and mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            found_groups.append((group_name, len(matches), group_info['losses']))

    # Generate fragments for each functional group
    seen_mz = set()  # Avoid duplicate m/z values

    for group_name, count, losses in found_groups:
        for loss_name, loss_mass in losses:
            fragment_mass = neutral_mass - loss_mass

            if fragment_mass < 50:  # Skip very small fragments
                continue

            # Calculate m/z based on charge
            if charge > 0:
                fragment_mz = (fragment_mass + charge * proton_mass) / abs(charge)
            else:
                fragment_mz = (fragment_mass - abs(charge) * proton_mass) / abs(charge)

            if fragment_mz <= 0:
                continue

            # Round to avoid floating point duplicates
            mz_key = round(fragment_mz, 4)
            if mz_key in seen_mz:
                continue
            seen_mz.add(mz_key)

            # Intensity scaling: more common groups get higher intensity
            relative_intensity = min(70.0, 30.0 + count * 10)

            fragments.append({
                'mz': fragment_mz,
                'formula': f'-{loss_name}',
                'loss': f'{group_name}:{loss_name}',
                'relative_intensity': relative_intensity
            })

    # Add ring-specific fragments for cyclic structures
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()

    if num_rings > 0:
        # Common ring cleavage patterns
        ring_losses = [
            ('C3H4', 40.031300),   # Cyclopropene equivalent
            ('C4H6', 54.046950),   # Butadiene (from ring opening)
            ('C5H6', 66.046950),   # Cyclopentadiene
            ('C6H6', 78.046950),   # Benzene
        ]

        for loss_name, loss_mass in ring_losses:
            fragment_mass = neutral_mass - loss_mass
            if fragment_mass < 50:
                continue

            if charge > 0:
                fragment_mz = (fragment_mass + charge * proton_mass) / abs(charge)
            else:
                fragment_mz = (fragment_mass - abs(charge) * proton_mass) / abs(charge)

            mz_key = round(fragment_mz, 4)
            if mz_key in seen_mz:
                continue
            seen_mz.add(mz_key)

            fragments.append({
                'mz': fragment_mz,
                'formula': f'-{loss_name}',
                'loss': f'ring:{loss_name}',
                'relative_intensity': 25.0
            })

    return fragments


def generate_fragment_library(formula: str, precursor_mz: float, charge: int = 1,
                              smiles: str = None) -> List[Dict]:
    """
    Generate theoretical fragment library from molecular formula and optionally SMILES.

    Args:
        formula: Molecular formula string
        precursor_mz: Precursor m/z
        charge: Ion charge state
        smiles: Optional SMILES string for structure-aware fragmentation

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

    # Add structure-aware fragments if SMILES is provided
    if smiles:
        structure_fragments = generate_structure_fragments(smiles, precursor_mz, charge)
        if structure_fragments:
            # Track existing m/z values to avoid duplicates
            existing_mz = {round(f['mz'], 4) for f in fragments}

            for frag in structure_fragments:
                mz_key = round(frag['mz'], 4)
                if mz_key not in existing_mz:
                    fragments.append(frag)
                    existing_mz.add(mz_key)

            print(f"  Added {len(structure_fragments)} structure-aware fragments from functional group analysis")

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

        # Generate theoretical fragments (with structure-aware fragmentation if SMILES available)
        smiles = ref_row.get('smiles', '') if not pd.isna(ref_row.get('smiles', '')) else ''
        theoretical_fragments = generate_fragment_library(
            ref_row['formula'],
            ref_row['precursor_mz'],
            ref_row['charge'],
            smiles=smiles
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
