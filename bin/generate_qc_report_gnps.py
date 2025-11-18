#!/usr/bin/env python3
"""
Generate QC Report for GNPS Reference Library Annotations

Creates HTML quality control report with visualizations for GNPS annotation workflow.

Author: MSBuddy Annotation Pipeline
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
import base64
from io import BytesIO
import sys

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 10


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate QC report for GNPS annotations"
    )
    parser.add_argument(
        '--annotations',
        type=str,
        required=True,
        help='Annotated peaks TSV file (with similarity metrics)'
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
        default='GNPS Library',
        help='Sample name for report header'
    )
    parser.add_argument(
        '--neighbors',
        type=str,
        default=None,
        help='Optional: Top neighbors TSV file for network visualization'
    )

    return parser.parse_args()


def fig_to_base64(fig) -> str:
    """
    Convert matplotlib figure to base64-encoded PNG.

    Args:
        fig: Matplotlib figure

    Returns:
        Base64-encoded PNG string
    """
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=100, bbox_inches='tight')
    buf.seek(0)
    img_str = base64.b64encode(buf.read()).decode('utf-8')
    buf.close()
    plt.close(fig)
    return img_str


def plot_annotation_quality_distribution(df: pd.DataFrame) -> str:
    """Plot annotation quality classification distribution."""
    if 'annotation_quality' not in df.columns:
        return ""

    fig, ax = plt.subplots(figsize=(8, 6))

    quality_counts = df['annotation_quality'].value_counts()

    # Color mapping
    colors = {
        'good': '#2ecc71',
        'uncertain': '#f39c12',
        'bad': '#e74c3c',
        'no_annotation': '#95a5a6'
    }

    # Get colors in order of quality_counts
    bar_colors = [colors.get(q, '#3498db') for q in quality_counts.index]

    bars = ax.bar(range(len(quality_counts)), quality_counts.values, color=bar_colors, alpha=0.7)
    ax.set_xticks(range(len(quality_counts)))
    ax.set_xticklabels(quality_counts.index, rotation=45, ha='right')
    ax.set_ylabel('Number of Spectra')
    ax.set_title('Annotation Quality Distribution', fontsize=12, fontweight='bold')

    # Add counts on bars
    for i, (bar, count) in enumerate(zip(bars, quality_counts.values)):
        height = bar.get_height()
        pct = count / len(df) * 100
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{count}\n({pct:.1f}%)',
                ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    return fig_to_base64(fig)


def plot_quality_score_distribution(df: pd.DataFrame) -> str:
    """Plot quality score distribution."""
    if 'quality_score' not in df.columns:
        return ""

    quality_scores = df['quality_score'].dropna()

    if len(quality_scores) == 0:
        return ""

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(quality_scores, bins=50, color='#3498db', alpha=0.7, edgecolor='black')
    ax.axvline(quality_scores.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {quality_scores.mean():.3f}')
    ax.axvline(quality_scores.median(), color='green', linestyle='--', linewidth=2, label=f'Median: {quality_scores.median():.3f}')
    ax.axvline(0.65, color='orange', linestyle=':', linewidth=2, label='Good Threshold (0.65, ~1% FDR)')
    ax.axvline(0.50, color='purple', linestyle=':', linewidth=2, label='Uncertain Threshold (0.50, ~5% FDR)')

    ax.set_xlabel('Quality Score')
    ax.set_ylabel('Number of Spectra')
    ax.set_title('Quality Score Distribution', fontsize=12, fontweight='bold')
    ax.legend()

    plt.tight_layout()
    return fig_to_base64(fig)


def plot_matched_peaks_distribution(df: pd.DataFrame) -> str:
    """Plot matched peaks percentage distribution."""
    if 'matched_peaks_pct' not in df.columns:
        return ""

    matched_pct = df['matched_peaks_pct'].dropna()

    if len(matched_pct) == 0:
        return ""

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(matched_pct, bins=50, color='#9b59b6', alpha=0.7, edgecolor='black')
    ax.axvline(matched_pct.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {matched_pct.mean():.1f}%')
    ax.axvline(matched_pct.median(), color='green', linestyle='--', linewidth=2, label=f'Median: {matched_pct.median():.1f}%')

    ax.set_xlabel('Matched Peaks (%)')
    ax.set_ylabel('Number of Spectra')
    ax.set_title('Matched Peaks Percentage Distribution', fontsize=12, fontweight='bold')
    ax.legend()

    plt.tight_layout()
    return fig_to_base64(fig)


def plot_matched_intensity_distribution(df: pd.DataFrame) -> str:
    """Plot matched intensity percentage distribution."""
    if 'matched_intensity_pct' not in df.columns:
        return ""

    intensity_pct = df['matched_intensity_pct'].dropna()

    if len(intensity_pct) == 0:
        return ""

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(intensity_pct, bins=50, color='#e67e22', alpha=0.7, edgecolor='black')
    ax.axvline(intensity_pct.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {intensity_pct.mean():.1f}%')
    ax.axvline(intensity_pct.median(), color='green', linestyle='--', linewidth=2, label=f'Median: {intensity_pct.median():.1f}%')

    ax.set_xlabel('Matched Intensity (%)')
    ax.set_ylabel('Number of Spectra')
    ax.set_title('Matched Intensity Percentage Distribution', fontsize=12, fontweight='bold')
    ax.legend()

    plt.tight_layout()
    return fig_to_base64(fig)


def plot_similarity_distribution(df: pd.DataFrame) -> str:
    """Plot cosine similarity score distribution."""
    if 'max_similarity_score' not in df.columns:
        return ""

    similarity_scores = df['max_similarity_score'].dropna()
    similarity_scores = similarity_scores[similarity_scores > 0]  # Remove zeros (no neighbors)

    if len(similarity_scores) == 0:
        return ""

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(similarity_scores, bins=50, color='#16a085', alpha=0.7, edgecolor='black')
    ax.axvline(similarity_scores.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {similarity_scores.mean():.3f}')
    ax.axvline(similarity_scores.median(), color='green', linestyle='--', linewidth=2, label=f'Median: {similarity_scores.median():.3f}')

    ax.set_xlabel('Maximum Cosine Similarity Score')
    ax.set_ylabel('Number of Spectra')
    ax.set_title('Spectral Similarity Distribution (Max Score per Spectrum)', fontsize=12, fontweight='bold')
    ax.legend()

    plt.tight_layout()
    return fig_to_base64(fig)


def plot_quality_vs_similarity(df: pd.DataFrame) -> str:
    """Plot quality score vs similarity score."""
    if 'quality_score' not in df.columns or 'max_similarity_score' not in df.columns:
        return ""

    # Filter out spectra with no neighbors
    plot_df = df[(df['quality_score'].notna()) & (df['max_similarity_score'] > 0)].copy()

    if len(plot_df) == 0:
        return ""

    fig, ax = plt.subplots(figsize=(10, 8))

    # Color by annotation quality
    if 'annotation_quality' in plot_df.columns:
        color_map = {
            'good': '#2ecc71',
            'uncertain': '#f39c12',
            'bad': '#e74c3c',
            'no_annotation': '#95a5a6'
        }
        colors = plot_df['annotation_quality'].map(color_map).fillna('#3498db')

        for quality in color_map.keys():
            mask = plot_df['annotation_quality'] == quality
            if mask.any():
                ax.scatter(plot_df.loc[mask, 'quality_score'],
                          plot_df.loc[mask, 'max_similarity_score'],
                          c=color_map[quality], label=quality, alpha=0.6, s=50)
    else:
        ax.scatter(plot_df['quality_score'], plot_df['max_similarity_score'],
                  alpha=0.6, s=50, c='#3498db')

    ax.set_xlabel('Quality Score (Peak Matching)')
    ax.set_ylabel('Max Similarity Score (Cosine)')
    ax.set_title('Quality Score vs. Spectral Similarity', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)

    if 'annotation_quality' in plot_df.columns:
        ax.legend(title='Annotation Quality')

    plt.tight_layout()
    return fig_to_base64(fig)


def plot_neighbors_distribution(df: pd.DataFrame) -> str:
    """Plot distribution of number of similar neighbors."""
    if 'num_similar_neighbors' not in df.columns:
        return ""

    neighbor_counts = df['num_similar_neighbors'].dropna()

    if len(neighbor_counts) == 0:
        return ""

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(neighbor_counts, bins=min(50, int(neighbor_counts.max())), color='#34495e', alpha=0.7, edgecolor='black')
    ax.axvline(neighbor_counts.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {neighbor_counts.mean():.1f}')
    ax.axvline(neighbor_counts.median(), color='green', linestyle='--', linewidth=2, label=f'Median: {neighbor_counts.median():.0f}')

    ax.set_xlabel('Number of Similar Neighbors')
    ax.set_ylabel('Number of Spectra')
    ax.set_title('Distribution of Similar Neighbors per Spectrum', fontsize=12, fontweight='bold')
    ax.legend()

    plt.tight_layout()
    return fig_to_base64(fig)


def plot_adduct_frequency(df: pd.DataFrame) -> str:
    """Plot adduct type frequency."""
    if 'adduct' not in df.columns:
        return ""

    adduct_counts = df['adduct'].value_counts().head(20)

    if len(adduct_counts) == 0:
        return ""

    fig, ax = plt.subplots(figsize=(10, 6))

    bars = ax.barh(range(len(adduct_counts)), adduct_counts.values, color='#2c3e50', alpha=0.7)
    ax.set_yticks(range(len(adduct_counts)))
    ax.set_yticklabels(adduct_counts.index)
    ax.set_xlabel('Number of Spectra')
    ax.set_title('Top 20 Adduct Types', fontsize=12, fontweight='bold')
    ax.invert_yaxis()

    # Add counts
    for i, (bar, count) in enumerate(zip(bars, adduct_counts.values)):
        width = bar.get_width()
        ax.text(width, bar.get_y() + bar.get_height()/2.,
                f' {count}',
                ha='left', va='center', fontsize=9)

    plt.tight_layout()
    return fig_to_base64(fig)


def generate_summary_table(df: pd.DataFrame) -> str:
    """Generate HTML summary statistics table."""
    total = len(df)

    table_html = """
    <table class="summary-table">
        <tr>
            <th colspan="2">Dataset Overview</th>
        </tr>
        <tr>
            <td>Total Spectra</td>
            <td>{total}</td>
        </tr>
    """.format(total=total)

    # Quality distribution
    if 'annotation_quality' in df.columns:
        quality_counts = df['annotation_quality'].value_counts()
        table_html += """
        <tr>
            <th colspan="2">Annotation Quality</th>
        </tr>
        """
        for quality in ['good', 'uncertain', 'bad', 'no_annotation']:
            count = quality_counts.get(quality, 0)
            pct = count / total * 100 if total > 0 else 0
            table_html += f"""
        <tr>
            <td>{quality.replace('_', ' ').title()}</td>
            <td>{count} ({pct:.1f}%)</td>
        </tr>
            """

    # Peak matching statistics
    if 'matched_peaks_pct' in df.columns:
        matched_pct = df['matched_peaks_pct'].dropna()
        if len(matched_pct) > 0:
            table_html += f"""
        <tr>
            <th colspan="2">Peak Matching</th>
        </tr>
        <tr>
            <td>Avg Matched Peaks</td>
            <td>{matched_pct.mean():.1f}%</td>
        </tr>
            """

    # Intensity matching statistics
    if 'matched_intensity_pct' in df.columns:
        intensity_pct = df['matched_intensity_pct'].dropna()
        if len(intensity_pct) > 0:
            table_html += f"""
        <tr>
            <td>Avg Matched Intensity</td>
            <td>{intensity_pct.mean():.1f}%</td>
        </tr>
            """

    # Quality score
    if 'quality_score' in df.columns:
        quality_score = df['quality_score'].dropna()
        if len(quality_score) > 0:
            table_html += f"""
        <tr>
            <td>Avg Quality Score</td>
            <td>{quality_score.mean():.3f}</td>
        </tr>
            """

    # Similarity statistics
    if 'max_similarity_score' in df.columns:
        sim_scores = df['max_similarity_score'].dropna()
        sim_scores_nonzero = sim_scores[sim_scores > 0]
        if len(sim_scores_nonzero) > 0:
            table_html += f"""
        <tr>
            <th colspan="2">Spectral Similarity</th>
        </tr>
        <tr>
            <td>Spectra with Neighbors</td>
            <td>{len(sim_scores_nonzero)} ({len(sim_scores_nonzero)/total*100:.1f}%)</td>
        </tr>
        <tr>
            <td>Avg Max Similarity</td>
            <td>{sim_scores_nonzero.mean():.3f}</td>
        </tr>
            """

    # Structural information
    if 'smiles' in df.columns:
        has_smiles = df['smiles'].notna() & (df['smiles'] != '')
        table_html += f"""
        <tr>
            <th colspan="2">Structural Information</th>
        </tr>
        <tr>
            <td>With SMILES</td>
            <td>{has_smiles.sum()} ({has_smiles.sum()/total*100:.1f}%)</td>
        </tr>
        """

    if 'formula' in df.columns:
        unique_formulas = df['formula'].dropna().nunique()
        table_html += f"""
        <tr>
            <td>Unique Formulas</td>
            <td>{unique_formulas}</td>
        </tr>
        """

    table_html += "</table>"
    return table_html


def generate_html_report(df: pd.DataFrame, sample_name: str) -> str:
    """Generate complete HTML report."""

    # Generate all plots
    plots = {
        'quality_dist': plot_annotation_quality_distribution(df),
        'quality_score': plot_quality_score_distribution(df),
        'matched_peaks': plot_matched_peaks_distribution(df),
        'matched_intensity': plot_matched_intensity_distribution(df),
        'similarity': plot_similarity_distribution(df),
        'quality_vs_sim': plot_quality_vs_similarity(df),
        'neighbors': plot_neighbors_distribution(df),
        'adducts': plot_adduct_frequency(df),
    }

    # Generate summary table
    summary_table = generate_summary_table(df)

    # Calculate annotation rate
    total = len(df)
    if 'annotation_quality' in df.columns:
        annotated = df[df['annotation_quality'] != 'no_annotation'].shape[0]
        annotation_rate = annotated / total * 100 if total > 0 else 0
    else:
        annotation_rate = 0

    # Build HTML
    html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GNPS Annotation QC Report - {sample_name}</title>
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
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}

        .header {{
            text-align: center;
            padding: 20px 0;
            border-bottom: 3px solid #3498db;
            margin-bottom: 30px;
        }}

        .header h1 {{
            color: #2c3e50;
            font-size: 2.5em;
            margin-bottom: 10px;
        }}

        .header .meta {{
            color: #7f8c8d;
            font-size: 1em;
        }}

        .annotation-rate {{
            text-align: center;
            padding: 30px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-radius: 8px;
            margin-bottom: 30px;
        }}

        .annotation-rate .value {{
            font-size: 4em;
            font-weight: bold;
            margin-bottom: 10px;
        }}

        .annotation-rate .label {{
            font-size: 1.2em;
            opacity: 0.9;
        }}

        .section {{
            margin: 40px 0;
        }}

        .section h2 {{
            color: #2c3e50;
            font-size: 1.8em;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid #ecf0f1;
        }}

        .summary-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}

        .summary-table th {{
            background-color: #3498db;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: 600;
        }}

        .summary-table td {{
            padding: 10px 12px;
            border-bottom: 1px solid #ecf0f1;
        }}

        .summary-table tr:nth-child(even) {{
            background-color: #f8f9fa;
        }}

        .plot {{
            margin: 30px 0;
            text-align: center;
        }}

        .plot h3 {{
            color: #34495e;
            margin-bottom: 15px;
            font-size: 1.3em;
        }}

        .plot img {{
            max-width: 100%;
            height: auto;
            border-radius: 4px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}

        .footer {{
            text-align: center;
            padding: 20px 0;
            margin-top: 40px;
            border-top: 2px solid #ecf0f1;
            color: #7f8c8d;
        }}

        @media print {{
            body {{
                background: white;
                padding: 0;
            }}
            .container {{
                box-shadow: none;
            }}
            .plot {{
                page-break-inside: avoid;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>GNPS Reference Library QC Report</h1>
            <div class="meta">
                <div>Sample: {sample_name}</div>
                <div>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
            </div>
        </div>

        <div class="annotation-rate">
            <div class="value">{annotation_rate:.1f}%</div>
            <div class="label">Annotation Success Rate</div>
        </div>

        <div class="section">
            <h2>Summary Statistics</h2>
            {summary_table}
        </div>

        <div class="section">
            <h2>Quality Control Visualizations</h2>
    """

    # Add plots
    if plots['quality_dist']:
        html += f"""
            <div class="plot">
                <h3>Annotation Quality Distribution</h3>
                <img src="data:image/png;base64,{plots['quality_dist']}" alt="Quality Distribution">
            </div>
        """

    if plots['quality_score']:
        html += f"""
            <div class="plot">
                <h3>Quality Score Distribution</h3>
                <img src="data:image/png;base64,{plots['quality_score']}" alt="Quality Score">
            </div>
        """

    if plots['matched_peaks']:
        html += f"""
            <div class="plot">
                <h3>Matched Peaks Percentage</h3>
                <img src="data:image/png;base64,{plots['matched_peaks']}" alt="Matched Peaks">
            </div>
        """

    if plots['matched_intensity']:
        html += f"""
            <div class="plot">
                <h3>Matched Intensity Percentage</h3>
                <img src="data:image/png;base64,{plots['matched_intensity']}" alt="Matched Intensity">
            </div>
        """

    if plots['similarity']:
        html += f"""
            <div class="plot">
                <h3>Spectral Similarity Distribution</h3>
                <img src="data:image/png;base64,{plots['similarity']}" alt="Similarity">
            </div>
        """

    if plots['quality_vs_sim']:
        html += f"""
            <div class="plot">
                <h3>Quality Score vs. Spectral Similarity</h3>
                <img src="data:image/png;base64,{plots['quality_vs_sim']}" alt="Quality vs Similarity">
            </div>
        """

    if plots['neighbors']:
        html += f"""
            <div class="plot">
                <h3>Similar Neighbors Distribution</h3>
                <img src="data:image/png;base64,{plots['neighbors']}" alt="Neighbors">
            </div>
        """

    if plots['adducts']:
        html += f"""
            <div class="plot">
                <h3>Adduct Type Frequency</h3>
                <img src="data:image/png;base64,{plots['adducts']}" alt="Adducts">
            </div>
        """

    html += """
        </div>

        <div class="footer">
            <p>GNPS Reference Library Annotation Pipeline v3.0</p>
            <p>Powered by MSBuddy-Annotation-Pipeline</p>
        </div>
    </div>
</body>
</html>
    """

    return html


def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 60)
    print("GNPS QC Report Generation")
    print("=" * 60)

    # Load annotations
    print(f"\nLoading annotations from: {args.annotations}")
    try:
        df = pd.read_csv(args.annotations, sep='\t')
        print(f"Loaded {len(df)} annotations")
    except Exception as e:
        print(f"Error loading annotations: {e}", file=sys.stderr)
        sys.exit(1)

    # Generate HTML report
    print("\nGenerating HTML report...")
    html_content = generate_html_report(df, args.sample_name)

    # Save report
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(html_content)

    print(f"Report saved to: {args.output}")
    print("=" * 60)


if __name__ == '__main__':
    main()
