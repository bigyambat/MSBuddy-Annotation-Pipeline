#!/usr/bin/env python3
"""
Parse GNPS Reference Library MGF Files

Extracts reference annotations (SMILES, InChI, molecular formula) from GNPS MGF files
and generates theoretical fragment predictions for peak annotation.

Author: MSBuddy Annotation Pipeline
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import sys
import re
from collections import defaultdict

try:
    from pyteomics import mgf
except ImportError:
    print("Error: pyteomics not installed. Install with: pip install pyteomics", file=sys.stderr)
    sys.exit(1)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse GNPS reference library MGF files and extract annotations"
    )
    parser.add_argument(
        '--mgf',
        type=str,
        required=True,
        help='Input MGF file with GNPS reference library data'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Output TSV file with reference annotations'
    )
    parser.add_argument(
        '--min_peaks',
        type=int,
        default=3,
        help='Minimum number of peaks required in spectrum (default: 3)'
    )

    return parser.parse_args()


def parse_molecular_formula(formula_str: str) -> Dict[str, int]:
    """
    Parse molecular formula string into element counts.

    Args:
        formula_str: Molecular formula (e.g., "C15H32O2")

    Returns:
        Dictionary mapping element symbols to counts

    Examples:
        >>> parse_molecular_formula("C15H32O2")
        {'C': 15, 'H': 32, 'O': 2}
    """
    if not formula_str or pd.isna(formula_str):
        return {}

    # Remove whitespace and common decorators
    formula_str = str(formula_str).strip()
    formula_str = re.sub(r'[\[\]]', '', formula_str)  # Remove brackets

    # Pattern: Element symbol (uppercase + optional lowercase) followed by optional number
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula_str)

    composition = {}
    for element, count in matches:
        if element:  # Skip empty matches
            count = int(count) if count else 1
            composition[element] = composition.get(element, 0) + count

    return composition


def calculate_formula_mass(composition: Dict[str, int]) -> float:
    """
    Calculate monoisotopic mass from molecular formula composition.

    Args:
        composition: Dictionary mapping element symbols to counts

    Returns:
        Monoisotopic mass in Daltons
    """
    # Monoisotopic masses (most abundant isotope)
    ATOMIC_MASSES = {
        'H': 1.007825032,
        'C': 12.0000000,
        'N': 14.003074005,
        'O': 15.994914620,
        'P': 30.973761998,
        'S': 31.972071174,
        'F': 18.998403163,
        'Cl': 34.968852682,
        'Br': 78.918337600,
        'I': 126.904473000,
        'Na': 22.989769282,
        'K': 38.963706487,
        'Ca': 39.962590863,
        'Mg': 23.985041697,
        'Fe': 55.934937475,
        'Se': 79.916521800,
        'Si': 27.976926535,
        'B': 11.009305167,
    }

    mass = 0.0
    for element, count in composition.items():
        if element in ATOMIC_MASSES:
            mass += ATOMIC_MASSES[element] * count
        else:
            print(f"Warning: Unknown element '{element}' in formula. Skipping.", file=sys.stderr)

    return mass


def generate_theoretical_fragments(composition: Dict[str, int],
                                   precursor_mz: float,
                                   charge: int = 1,
                                   max_fragments: int = 20) -> List[Dict]:
    """
    Generate theoretical fragment m/z values from molecular formula.

    This is a simplified approach generating common neutral losses and fragments.
    For GNPS reference library annotation, we'll primarily rely on matching
    observed peaks to theoretical possibilities.

    Args:
        composition: Molecular formula as element counts
        precursor_mz: Precursor m/z value
        charge: Ion charge state
        max_fragments: Maximum number of theoretical fragments to generate

    Returns:
        List of theoretical fragments with m/z and formula
    """
    if not composition:
        return []

    fragments = []
    proton_mass = 1.007276467
    electron_mass = 0.000548579909

    # Common neutral losses (in Da)
    COMMON_LOSSES = {
        'H2O': 18.010565,      # Water loss
        'NH3': 17.026549,      # Ammonia loss
        'CO': 27.994915,       # CO loss
        'CO2': 43.989829,      # CO2 loss
        'CH3': 15.023475,      # Methyl loss
        'C2H4': 28.031300,     # Ethylene loss
        'H2': 2.015650,        # Hydrogen loss
    }

    # Calculate neutral mass
    if charge > 0:  # Positive mode
        neutral_mass = (precursor_mz * abs(charge)) - (charge * proton_mass)
    else:  # Negative mode
        neutral_mass = (precursor_mz * abs(charge)) + (abs(charge) * (proton_mass - 2 * electron_mass))

    # Add precursor ion
    fragments.append({
        'mz': precursor_mz,
        'formula': 'M',
        'loss': 'precursor',
        'intensity': 100.0  # Placeholder
    })

    # Generate common neutral loss fragments
    for loss_name, loss_mass in COMMON_LOSSES.items():
        # Check if loss is chemically possible
        if loss_name == 'H2O' and composition.get('H', 0) >= 2 and composition.get('O', 0) >= 1:
            fragment_mass = neutral_mass - loss_mass
        elif loss_name == 'NH3' and composition.get('N', 0) >= 1 and composition.get('H', 0) >= 3:
            fragment_mass = neutral_mass - loss_mass
        elif loss_name == 'CO2' and composition.get('C', 0) >= 1 and composition.get('O', 0) >= 2:
            fragment_mass = neutral_mass - loss_mass
        elif loss_name == 'CO' and composition.get('C', 0) >= 1 and composition.get('O', 0) >= 1:
            fragment_mass = neutral_mass - loss_mass
        else:
            continue

        # Convert to m/z (assume same charge as precursor)
        if charge > 0:
            fragment_mz = (fragment_mass + charge * proton_mass) / abs(charge)
        else:
            fragment_mz = (fragment_mass - abs(charge) * (proton_mass - 2 * electron_mass)) / abs(charge)

        if fragment_mz > 0:  # Valid m/z
            fragments.append({
                'mz': fragment_mz,
                'formula': f'M-{loss_name}',
                'loss': loss_name,
                'intensity': 50.0  # Placeholder
            })

    return fragments[:max_fragments]


def parse_gnps_mgf(mgf_path: str, min_peaks: int = 3) -> pd.DataFrame:
    """
    Parse GNPS reference library MGF file and extract annotations.

    Args:
        mgf_path: Path to MGF file
        min_peaks: Minimum number of peaks required

    Returns:
        DataFrame with columns:
            - spectrum_id: Unique spectrum identifier
            - smiles: SMILES string
            - inchi: InChI string
            - inchikey: InChIKey
            - formula: Molecular formula
            - precursor_mz: Precursor m/z
            - charge: Ion charge
            - adduct: Adduct type
            - num_peaks: Number of peaks in spectrum
            - formula_mass: Calculated mass from formula
            - collision_energy: Collision energy (if available)
            - retention_time: Retention time (if available)
            - compound_name: Compound name (if available)
            - library_id: GNPS library ID
    """
    records = []
    skipped = 0

    print(f"Reading MGF file: {mgf_path}")

    try:
        with mgf.read(mgf_path, use_index=False) as reader:
            for idx, spectrum in enumerate(reader):
                # Extract spectrum parameters
                params = spectrum.get('params', {})

                # Extract spectrum ID (priority order)
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

                # Skip if missing required data
                if mz_array is None or intensity_array is None:
                    skipped += 1
                    continue

                if len(mz_array) < min_peaks:
                    skipped += 1
                    continue

                # Extract GNPS-specific fields
                smiles = params.get('smiles', '')
                inchi = params.get('inchi', '')
                inchikey = params.get('inchikey', '')
                formula = params.get('formula', '') or params.get('molecularformula', '')
                compound_name = params.get('name', '') or params.get('compound_name', '')
                library_id = params.get('spectrumid', '')

                # Extract precursor information
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

                # Extract adduct
                adduct = params.get('adduct', '') or params.get('ionmode', '')

                # Extract optional metadata
                collision_energy = params.get('collision_energy', '') or params.get('collisionenergy', '')
                retention_time = params.get('rtinseconds', '') or params.get('retention_time', '')

                # Parse molecular formula and calculate mass
                composition = parse_molecular_formula(formula)
                formula_mass = calculate_formula_mass(composition) if composition else 0.0

                # Create record
                record = {
                    'spectrum_id': spec_id,
                    'smiles': smiles,
                    'inchi': inchi,
                    'inchikey': inchikey,
                    'formula': formula,
                    'precursor_mz': precursor_mz,
                    'charge': charge,
                    'adduct': adduct,
                    'num_peaks': len(mz_array),
                    'formula_mass': formula_mass,
                    'collision_energy': collision_energy,
                    'retention_time': retention_time,
                    'compound_name': compound_name,
                    'library_id': library_id,
                }

                records.append(record)

                # Progress update
                if (idx + 1) % 1000 == 0:
                    print(f"Processed {idx + 1} spectra...")

    except Exception as e:
        print(f"Error reading MGF file: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Parsed {len(records)} spectra (skipped {skipped} with insufficient peaks)")

    # Create DataFrame
    df = pd.DataFrame(records)

    return df


def validate_annotations(df: pd.DataFrame) -> pd.DataFrame:
    """
    Validate and report on annotation quality.

    Args:
        df: DataFrame with reference annotations

    Returns:
        DataFrame with validation flags added
    """
    print("\n=== Annotation Validation Summary ===")

    total = len(df)

    # Check for required fields
    has_smiles = df['smiles'].notna() & (df['smiles'] != '')
    has_inchi = df['inchi'].notna() & (df['inchi'] != '')
    has_formula = df['formula'].notna() & (df['formula'] != '')
    has_name = df['compound_name'].notna() & (df['compound_name'] != '')

    print(f"Total spectra: {total}")
    print(f"With SMILES: {has_smiles.sum()} ({has_smiles.sum()/total*100:.1f}%)")
    print(f"With InChI: {has_inchi.sum()} ({has_inchi.sum()/total*100:.1f}%)")
    print(f"With Formula: {has_formula.sum()} ({has_formula.sum()/total*100:.1f}%)")
    print(f"With Name: {has_name.sum()} ({has_name.sum()/total*100:.1f}%)")

    # Add validation flags
    df['has_structure'] = has_smiles | has_inchi
    df['has_formula'] = has_formula
    df['annotation_complete'] = has_smiles & has_formula

    # Check mass accuracy (if formula available)
    if has_formula.any():
        df['mass_error_ppm'] = np.where(
            (df['formula_mass'] > 0) & (df['precursor_mz'] > 0),
            abs(df['formula_mass'] - df['precursor_mz']) / df['precursor_mz'] * 1e6,
            np.nan
        )

        valid_errors = df['mass_error_ppm'].dropna()
        if len(valid_errors) > 0:
            print(f"\nMass Accuracy:")
            print(f"  Mean error: {valid_errors.mean():.2f} ppm")
            print(f"  Median error: {valid_errors.median():.2f} ppm")
            print(f"  Std dev: {valid_errors.std():.2f} ppm")

    # Report annotation completeness
    complete = df['annotation_complete'].sum()
    print(f"\nComplete annotations (SMILES + Formula): {complete} ({complete/total*100:.1f}%)")

    return df


def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 60)
    print("GNPS Reference Library Parser")
    print("=" * 60)

    # Parse MGF file
    df = parse_gnps_mgf(args.mgf, min_peaks=args.min_peaks)

    if len(df) == 0:
        print("Error: No valid spectra found in MGF file", file=sys.stderr)
        sys.exit(1)

    # Validate annotations
    df = validate_annotations(df)

    # Save to TSV
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(args.output, sep='\t', index=False)
    print(f"\nReference annotations saved to: {args.output}")
    print(f"Total records: {len(df)}")
    print("=" * 60)


if __name__ == '__main__':
    main()
