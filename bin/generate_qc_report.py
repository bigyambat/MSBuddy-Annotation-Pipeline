#!/usr/bin/env python3

"""
Mass Spectrum Annotation QC Report Generator

Author: Bigy Ambat
Version: 2.0
Description: Generates comprehensive HTML QC reports for MSBuddy annotation results
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from pyteomics import mgf
import base64
from io import BytesIO
from datetime import datetime

# Set style for plots
sns.set_theme(style="whitegrid")
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 10


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Generate QC report from MSBuddy annotation results'
    )
    parser.add_argument(
        '--mgf',
        type=str,
        required=True,
        help='Input MGF file'
    )
    parser.add_argument(
        '--tsv',
        type=str,
        required=True,
        help='MSBuddy annotation TSV file'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Output HTML report file'
    )
    parser.add_argument(
        '--sample_name',
        type=str,
        required=True,
        help='Sample name for the report'
    )
    return parser.parse_args()


def count_spectra_from_mgf(mgf_file):
    """
    Count total number of spectra in MGF file using pyteomics

    Args:
        mgf_file (str): Path to MGF file

    Returns:
        int: Number of spectra in the file
    """
    try:
        with mgf.read(mgf_file) as reader:
            count = sum(1 for _ in reader)
        return count
    except Exception as e:
        print(f"Error reading MGF file: {e}", file=sys.stderr)
        sys.exit(1)


def load_annotation_results(tsv_file):
    """
    Load MSBuddy annotation results from TSV file

    Args:
        tsv_file (str): Path to TSV annotation file

    Returns:
        pd.DataFrame: Annotation results
    """
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        return df
    except Exception as e:
        print(f"Error reading annotation TSV file: {e}", file=sys.stderr)
        sys.exit(1)


def calculate_annotation_rate(total_spectra, df):
    """
    Calculate annotation rate

    Args:
        total_spectra (int): Total number of spectra
        df (pd.DataFrame): Annotation results dataframe

    Returns:
        tuple: (annotated_count, annotation_rate)
    """
    # Count non-null annotations (assuming annotations have formula or structure info)
    # Try multiple possible column names that MSBuddy might use
    possible_formula_cols = [
        'formula_rank_1', 'formula', 'molecular_formula', 'mf',
        'predicted_formula', 'top_formula', 'formula_1'
    ]

    annotated = df
    found_col = None

    for col in possible_formula_cols:
        if col in df.columns:
            # Count rows where the column has a non-empty, non-null value
            annotated = df[df[col].notna() & (df[col] != '') & (df[col] != 'None')]
            found_col = col
            break

    # If no formula column found, check for other indicators of successful annotation
    if found_col is None:
        # Check for score column - if present and > 0, likely annotated
        score_cols = ['score', 'msbuddy_score', 'confidence', 'final_score']
        for col in score_cols:
            if col in df.columns:
                annotated = df[df[col].notna() & (df[col] > 0)]
                found_col = col
                break

    # If still no indicator found, check if annotation_quality column exists (added by analyze_peak_explanation.py)
    if found_col is None and 'annotation_quality' in df.columns:
        # Count anything that's not "no_annotation"
        annotated = df[df['annotation_quality'] != 'no_annotation']
        found_col = 'annotation_quality'

    # Last resort: if we found nothing, count rows, but this is likely wrong
    if found_col is None:
        print(f"WARNING: Could not find formula or score columns. Available columns: {list(df.columns)}")
        print(f"WARNING: Assuming all {len(df)} rows are annotated (may be incorrect)")
        annotated = df

    annotated_count = len(annotated)
    annotation_rate = (annotated_count / total_spectra * 100) if total_spectra > 0 else 0

    print(f"Used column '{found_col}' to determine annotation rate")
    return annotated_count, annotation_rate


def fig_to_base64(fig):
    """
    Convert matplotlib figure to base64 encoded string for embedding in HTML

    Args:
        fig: matplotlib figure object

    Returns:
        str: Base64 encoded image string
    """
    buffer = BytesIO()
    fig.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
    buffer.seek(0)
    image_base64 = base64.b64encode(buffer.read()).decode()
    plt.close(fig)
    return image_base64


def plot_score_distribution(df):
    """
    Generate histogram of MSBuddy annotation scores

    Args:
        df (pd.DataFrame): Annotation results

    Returns:
        str: Base64 encoded plot
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Assuming MSBuddy has a 'score' column - adjust as needed
    score_col = 'score' if 'score' in df.columns else 'msbuddy_score'

    if score_col in df.columns:
        scores = df[score_col].dropna()
        ax.hist(scores, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
        ax.set_xlabel('MSBuddy Score', fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)
        ax.set_title('Distribution of MSBuddy Annotation Scores', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)

        # Add statistics
        mean_score = scores.mean()
        median_score = scores.median()
        ax.axvline(mean_score, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_score:.2f}')
        ax.axvline(median_score, color='green', linestyle='--', linewidth=2, label=f'Median: {median_score:.2f}')
        ax.legend()
    else:
        ax.text(0.5, 0.5, 'Score data not available',
                ha='center', va='center', fontsize=14, transform=ax.transAxes)

    return fig_to_base64(fig)


def plot_mass_error_distribution(df):
    """
    Generate histogram of precursor mass error in ppm

    Args:
        df (pd.DataFrame): Annotation results

    Returns:
        str: Base64 encoded plot
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Common column names for mass error
    mass_error_cols = ['mass_error_ppm', 'ppm_error', 'mass_error', 'delta_ppm']
    mass_error_col = None

    for col in mass_error_cols:
        if col in df.columns:
            mass_error_col = col
            break

    if mass_error_col:
        mass_errors = df[mass_error_col].dropna()
        ax.hist(mass_errors, bins=50, edgecolor='black', alpha=0.7, color='coral')
        ax.set_xlabel('Mass Error (ppm)', fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)
        ax.set_title('Distribution of Precursor Mass Error', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)

        # Add statistics
        mean_error = mass_errors.mean()
        median_error = mass_errors.median()
        ax.axvline(mean_error, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_error:.2f} ppm')
        ax.axvline(median_error, color='green', linestyle='--', linewidth=2, label=f'Median: {median_error:.2f} ppm')
        ax.legend()
    else:
        ax.text(0.5, 0.5, 'Mass error data not available',
                ha='center', va='center', fontsize=14, transform=ax.transAxes)

    return fig_to_base64(fig)


def plot_adduct_frequency(df):
    """
    Generate bar chart of identified adduct frequency

    Args:
        df (pd.DataFrame): Annotation results

    Returns:
        str: Base64 encoded plot
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    adduct_col = 'adduct' if 'adduct' in df.columns else 'adduct_type'

    if adduct_col in df.columns:
        adduct_counts = df[adduct_col].value_counts().head(20)  # Top 20 adducts

        adduct_counts.plot(kind='bar', ax=ax, color='seagreen', edgecolor='black', alpha=0.7)
        ax.set_xlabel('Adduct Type', fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)
        ax.set_title('Frequency of Identified Adducts', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        plt.xticks(rotation=45, ha='right')

        # Add value labels on bars
        for i, v in enumerate(adduct_counts.values):
            ax.text(i, v + max(adduct_counts.values) * 0.01, str(v),
                   ha='center', va='bottom', fontsize=9)
    else:
        ax.text(0.5, 0.5, 'Adduct data not available',
                ha='center', va='center', fontsize=14, transform=ax.transAxes)

    plt.tight_layout()
    return fig_to_base64(fig)


def plot_annotation_quality_distribution(df):
    """
    Plot distribution of annotation quality classifications (good/bad/uncertain)

    Args:
        df (pd.DataFrame): Annotation results with 'annotation_quality' column

    Returns:
        str: Base64 encoded PNG image
    """
    if 'annotation_quality' not in df.columns:
        return None

    plt.figure(figsize=(8, 6))

    quality_counts = df['annotation_quality'].value_counts()

    # Define colors for each category
    colors = {
        'good': '#4CAF50',
        'uncertain': '#FFC107',
        'bad': '#F44336',
        'no_annotation': '#9E9E9E'
    }
    plot_colors = [colors.get(cat, '#9E9E9E') for cat in quality_counts.index]

    bars = plt.bar(quality_counts.index, quality_counts.values, color=plot_colors, edgecolor='black', linewidth=1.2)

    plt.xlabel('Annotation Quality', fontsize=12, fontweight='bold')
    plt.ylabel('Count', fontsize=12, fontweight='bold')
    plt.title('Annotation Quality Classification', fontsize=14, fontweight='bold', pad=15)

    # Add count labels on bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}\n({height/len(df)*100:.1f}%)',
                ha='center', va='bottom', fontsize=10, fontweight='bold')

    plt.xticks(rotation=0, fontsize=11)
    plt.tight_layout()

    # Convert to base64
    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
    buffer.seek(0)
    image_base64 = base64.b64encode(buffer.read()).decode()
    plt.close()

    return f"data:image/png;base64,{image_base64}"


def plot_explained_peaks_distribution(df):
    """
    Plot distribution of explained peaks percentage

    Args:
        df (pd.DataFrame): Annotation results with 'explained_peaks_pct' column

    Returns:
        str: Base64 encoded PNG image
    """
    if 'explained_peaks_pct' not in df.columns:
        return None

    # Filter out zero values for better visualization
    data = df[df['explained_peaks_pct'] > 0]['explained_peaks_pct']

    if len(data) == 0:
        return None

    plt.figure(figsize=(10, 6))

    plt.hist(data, bins=30, color='#2196F3', edgecolor='black', alpha=0.7)

    mean_val = data.mean()
    median_val = data.median()

    plt.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.1f}%')
    plt.axvline(median_val, color='green', linestyle='--', linewidth=2, label=f'Median: {median_val:.1f}%')

    plt.xlabel('Explained Peaks (%)', fontsize=12, fontweight='bold')
    plt.ylabel('Frequency', fontsize=12, fontweight='bold')
    plt.title('Distribution of Explained Peaks Percentage', fontsize=14, fontweight='bold', pad=15)
    plt.legend(fontsize=11)
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()

    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
    buffer.seek(0)
    image_base64 = base64.b64encode(buffer.read()).decode()
    plt.close()

    return f"data:image/png;base64,{image_base64}"


def plot_explained_intensity_distribution(df):
    """
    Plot distribution of explained intensity percentage

    Args:
        df (pd.DataFrame): Annotation results with 'explained_intensity_pct' column

    Returns:
        str: Base64 encoded PNG image
    """
    if 'explained_intensity_pct' not in df.columns:
        return None

    # Filter out zero values
    data = df[df['explained_intensity_pct'] > 0]['explained_intensity_pct']

    if len(data) == 0:
        return None

    plt.figure(figsize=(10, 6))

    plt.hist(data, bins=30, color='#9C27B0', edgecolor='black', alpha=0.7)

    mean_val = data.mean()
    median_val = data.median()

    plt.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.1f}%')
    plt.axvline(median_val, color='green', linestyle='--', linewidth=2, label=f'Median: {median_val:.1f}%')

    plt.xlabel('Explained Intensity (%)', fontsize=12, fontweight='bold')
    plt.ylabel('Frequency', fontsize=12, fontweight='bold')
    plt.title('Distribution of Explained Intensity Percentage', fontsize=14, fontweight='bold', pad=15)
    plt.legend(fontsize=11)
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()

    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
    buffer.seek(0)
    image_base64 = base64.b64encode(buffer.read()).decode()
    plt.close()

    return f"data:image/png;base64,{image_base64}"


def generate_summary_table(total_spectra, annotated_count, annotation_rate, df):
    """
    Generate HTML summary statistics table

    Args:
        total_spectra (int): Total spectra count
        annotated_count (int): Number of annotated spectra
        annotation_rate (float): Annotation rate percentage
        df (pd.DataFrame): Annotation results

    Returns:
        str: HTML table string
    """
    # Calculate additional statistics
    if 'formula_rank_1' in df.columns:
        unique_formulas = df['formula_rank_1'].nunique()
    elif 'formula' in df.columns:
        unique_formulas = df['formula'].nunique()
    else:
        unique_formulas = 'N/A'

    score_col = 'score' if 'score' in df.columns else 'msbuddy_score'
    avg_score = f"{df[score_col].mean():.2f}" if score_col in df.columns else 'N/A'

    mass_error_cols = ['mass_error_ppm', 'ppm_error', 'mass_error', 'delta_ppm']
    mass_error_col = next((col for col in mass_error_cols if col in df.columns), None)
    avg_mass_error = f"{df[mass_error_col].mean():.2f}" if mass_error_col else 'N/A'

    adduct_col = 'adduct' if 'adduct' in df.columns else 'adduct_type'
    unique_adducts = df[adduct_col].nunique() if adduct_col in df.columns else 'N/A'

    # Add quality classification statistics if available
    quality_stats = ""
    if 'annotation_quality' in df.columns:
        quality_counts = df['annotation_quality'].value_counts()
        good_count = quality_counts.get('good', 0)
        bad_count = quality_counts.get('bad', 0)
        uncertain_count = quality_counts.get('uncertain', 0)
        good_pct = (good_count / len(df) * 100) if len(df) > 0 else 0

        quality_stats = f"""
        <tr>
            <td><strong>Good Quality Annotations</strong></td>
            <td><span style="color: #4CAF50; font-weight: bold;">{good_count:,} ({good_pct:.1f}%)</span></td>
        </tr>
        <tr>
            <td><strong>Uncertain Quality Annotations</strong></td>
            <td><span style="color: #FFC107; font-weight: bold;">{uncertain_count:,}</span></td>
        </tr>
        <tr>
            <td><strong>Poor Quality Annotations</strong></td>
            <td><span style="color: #F44336; font-weight: bold;">{bad_count:,}</span></td>
        </tr>
        """

    # Add peak explanation statistics if available
    peak_stats = ""
    if 'explained_peaks_pct' in df.columns:
        valid_data = df[df['explained_peaks_pct'] > 0]['explained_peaks_pct']
        if len(valid_data) > 0:
            avg_explained_peaks = f"{valid_data.mean():.1f}%"
            median_explained_peaks = f"{valid_data.median():.1f}%"
        else:
            avg_explained_peaks = "N/A"
            median_explained_peaks = "N/A"

        peak_stats += f"""
        <tr>
            <td>Average Explained Peaks</td>
            <td>{avg_explained_peaks}</td>
        </tr>
        <tr>
            <td>Median Explained Peaks</td>
            <td>{median_explained_peaks}</td>
        </tr>
        """

    if 'explained_intensity_pct' in df.columns:
        valid_data = df[df['explained_intensity_pct'] > 0]['explained_intensity_pct']
        if len(valid_data) > 0:
            avg_explained_intensity = f"{valid_data.mean():.1f}%"
            median_explained_intensity = f"{valid_data.median():.1f}%"
        else:
            avg_explained_intensity = "N/A"
            median_explained_intensity = "N/A"

        peak_stats += f"""
        <tr>
            <td>Average Explained Intensity</td>
            <td>{avg_explained_intensity}</td>
        </tr>
        <tr>
            <td>Median Explained Intensity</td>
            <td>{median_explained_intensity}</td>
        </tr>
        """

    html = f"""
    <table class="summary-table">
        <tr>
            <th>Metric</th>
            <th>Value</th>
        </tr>
        <tr>
            <td>Total Spectra</td>
            <td>{total_spectra:,}</td>
        </tr>
        <tr>
            <td>Annotated Spectra</td>
            <td>{annotated_count:,}</td>
        </tr>
        <tr>
            <td>Annotation Rate</td>
            <td><strong>{annotation_rate:.2f}%</strong></td>
        </tr>
        {quality_stats}
        <tr>
            <td>Unique Molecular Formulas</td>
            <td>{unique_formulas}</td>
        </tr>
        <tr>
            <td>Average Score</td>
            <td>{avg_score}</td>
        </tr>
        <tr>
            <td>Average Mass Error (ppm)</td>
            <td>{avg_mass_error}</td>
        </tr>
        {peak_stats}
        <tr>
            <td>Unique Adducts Identified</td>
            <td>{unique_adducts}</td>
        </tr>
    </table>
    """
    return html


def generate_html_report(sample_name, total_spectra, annotated_count, annotation_rate,
                        df, score_plot, mass_error_plot, adduct_plot, quality_plot=None,
                        explained_peaks_plot=None, explained_intensity_plot=None):
    """
    Generate complete HTML report

    Args:
        sample_name (str): Sample name
        total_spectra (int): Total spectra count
        annotated_count (int): Annotated spectra count
        annotation_rate (float): Annotation rate percentage
        df (pd.DataFrame): Annotation results
        score_plot (str): Base64 encoded score plot
        mass_error_plot (str): Base64 encoded mass error plot
        adduct_plot (str): Base64 encoded adduct plot
        quality_plot (str, optional): Base64 encoded quality classification plot
        explained_peaks_plot (str, optional): Base64 encoded explained peaks plot
        explained_intensity_plot (str, optional): Base64 encoded explained intensity plot

    Returns:
        str: Complete HTML report
    """
    summary_table = generate_summary_table(total_spectra, annotated_count, annotation_rate, df)
    generation_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>QC Report - {sample_name}</title>
        <style>
            * {{
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }}

            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                line-height: 1.6;
                color: #333;
                background-color: #f4f4f4;
                padding: 20px;
            }}

            .container {{
                max-width: 1400px;
                margin: 0 auto;
                background-color: white;
                padding: 40px;
                box-shadow: 0 0 20px rgba(0,0,0,0.1);
                border-radius: 8px;
            }}

            .header {{
                text-align: center;
                padding-bottom: 30px;
                border-bottom: 3px solid #4CAF50;
                margin-bottom: 30px;
            }}

            .header h1 {{
                color: #2c3e50;
                font-size: 2.5em;
                margin-bottom: 10px;
            }}

            .header .sample-name {{
                font-size: 1.5em;
                color: #555;
                font-weight: 300;
            }}

            .header .date {{
                color: #777;
                font-size: 0.9em;
                margin-top: 10px;
            }}

            .section {{
                margin-bottom: 40px;
            }}

            .section h2 {{
                color: #2c3e50;
                font-size: 1.8em;
                margin-bottom: 20px;
                padding-bottom: 10px;
                border-bottom: 2px solid #e0e0e0;
            }}

            .summary-table {{
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
                font-size: 1.1em;
                box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            }}

            .summary-table th {{
                background-color: #4CAF50;
                color: white;
                text-align: left;
                padding: 15px;
                font-weight: 600;
            }}

            .summary-table td {{
                padding: 12px 15px;
                border-bottom: 1px solid #ddd;
            }}

            .summary-table tr:nth-child(even) {{
                background-color: #f9f9f9;
            }}

            .summary-table tr:hover {{
                background-color: #f5f5f5;
            }}

            .plot {{
                text-align: center;
                margin: 30px 0;
                padding: 20px;
                background-color: #fafafa;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.05);
            }}

            .plot img {{
                max-width: 100%;
                height: auto;
                border-radius: 4px;
            }}

            .plot h3 {{
                color: #2c3e50;
                margin-bottom: 15px;
                font-size: 1.3em;
            }}

            .annotation-rate {{
                text-align: center;
                padding: 30px;
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                border-radius: 8px;
                margin: 20px 0;
                box-shadow: 0 4px 12px rgba(0,0,0,0.15);
            }}

            .annotation-rate .rate {{
                font-size: 4em;
                font-weight: bold;
                margin: 10px 0;
            }}

            .annotation-rate .label {{
                font-size: 1.3em;
                opacity: 0.95;
            }}

            .footer {{
                text-align: center;
                margin-top: 50px;
                padding-top: 20px;
                border-top: 2px solid #e0e0e0;
                color: #777;
                font-size: 0.9em;
            }}

            .grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
                gap: 20px;
                margin: 20px 0;
            }}

            @media print {{
                body {{
                    background-color: white;
                }}
                .container {{
                    box-shadow: none;
                }}
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>Mass Spectrum QC Report</h1>
                <div class="sample-name">{sample_name}</div>
                <div class="date">Generated: {generation_date}</div>
            </div>

            <div class="annotation-rate">
                <div class="label">Overall Annotation Rate</div>
                <div class="rate">{annotation_rate:.2f}%</div>
                <div class="label">{annotated_count:,} of {total_spectra:,} spectra annotated</div>
            </div>

            <div class="section">
                <h2>Summary Statistics</h2>
                {summary_table}
            </div>

            <div class="section">
                <h2>Quality Control Visualizations</h2>

                <div class="plot">
                    <h3>MSBuddy Score Distribution</h3>
                    <img src="data:image/png;base64,{score_plot}" alt="Score Distribution">
                </div>

                <div class="plot">
                    <h3>Precursor Mass Error Distribution</h3>
                    <img src="data:image/png;base64,{mass_error_plot}" alt="Mass Error Distribution">
                </div>

                <div class="plot">
                    <h3>Adduct Frequency</h3>
                    <img src="data:image/png;base64,{adduct_plot}" alt="Adduct Frequency">
                </div>
                """

    # Add quality classification plot if available
    if quality_plot:
        html += f"""
                <div class="plot">
                    <h3>Annotation Quality Classification</h3>
                    <img src="{quality_plot}" alt="Quality Classification">
                </div>
        """

    # Add explained peaks plot if available
    if explained_peaks_plot:
        html += f"""
                <div class="plot">
                    <h3>Explained Peaks Distribution</h3>
                    <img src="{explained_peaks_plot}" alt="Explained Peaks Distribution">
                </div>
        """

    # Add explained intensity plot if available
    if explained_intensity_plot:
        html += f"""
                <div class="plot">
                    <h3>Explained Intensity Distribution</h3>
                    <img src="{explained_intensity_plot}" alt="Explained Intensity Distribution">
                </div>
        """

    html += """
            </div>

            <div class="footer">
                <p>MS Annotation & QC Pipeline v2.1</p>
                <p>Powered by MSBuddy, Nextflow, and Python</p>
            </div>
        </div>
    </body>
    </html>
    """
    return html


def main():
    """Main function"""
    args = parse_args()

    print(f"Starting QC report generation for {args.sample_name}")
    print(f"MGF file: {args.mgf}")
    print(f"TSV file: {args.tsv}")

    # Count total spectra from MGF
    print("Counting spectra in MGF file...")
    total_spectra = count_spectra_from_mgf(args.mgf)
    print(f"Total spectra: {total_spectra}")

    # Load annotation results
    print("Loading annotation results...")
    df = load_annotation_results(args.tsv)
    print(f"Loaded {len(df)} annotation records")

    # Calculate annotation rate
    annotated_count, annotation_rate = calculate_annotation_rate(total_spectra, df)
    print(f"Annotation rate: {annotation_rate:.2f}%")

    # Generate plots
    print("Generating visualizations...")
    score_plot = plot_score_distribution(df)
    mass_error_plot = plot_mass_error_distribution(df)
    adduct_plot = plot_adduct_frequency(df)

    # Generate enhanced plots if columns exist
    quality_plot = None
    explained_peaks_plot = None
    explained_intensity_plot = None

    if 'annotation_quality' in df.columns:
        print("Generating annotation quality visualization...")
        quality_plot = plot_annotation_quality_distribution(df)

    if 'explained_peaks_pct' in df.columns:
        print("Generating explained peaks visualization...")
        explained_peaks_plot = plot_explained_peaks_distribution(df)

    if 'explained_intensity_pct' in df.columns:
        print("Generating explained intensity visualization...")
        explained_intensity_plot = plot_explained_intensity_distribution(df)

    # Generate HTML report
    print("Generating HTML report...")
    html_report = generate_html_report(
        args.sample_name,
        total_spectra,
        annotated_count,
        annotation_rate,
        df,
        score_plot,
        mass_error_plot,
        adduct_plot,
        quality_plot,
        explained_peaks_plot,
        explained_intensity_plot
    )

    # Write report to file
    with open(args.output, 'w') as f:
        f.write(html_report)

    print(f"QC report successfully generated: {args.output}")


if __name__ == '__main__':
    main()
