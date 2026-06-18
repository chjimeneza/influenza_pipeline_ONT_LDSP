#!/usr/bin/env python3
"""
Generate plots and analysis tables from influenza pipeline results.
Run this script after the main pipeline completes.
"""

import os
import subprocess
import re
from pathlib import Path
from collections import Counter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Configuration
REMAP_DIR = "remap_results"
REFS_DIR = "refs_by_barcode"
MULTIFASTA_DIR = "multi_fasta_segments"
REFERENCE_FASTA = "consensus_irma.fasta"
OUTPUT_PREFIX = "analysis"

# Segment list (match the pipeline)
SEGMENT_LIST = [
    'A_MP', 'A_NP', 'A_NS', 'A_PA', 'A_PB1', 'A_PB2',
    'A_HA_H1', 'A_HA_H10', 'A_HA_H11', 'A_HA_H12', 'A_HA_H13',
    'A_HA_H14', 'A_HA_H15', 'A_HA_H16', 'A_HA_H2', 'A_HA_H3',
    'A_HA_H4', 'A_HA_H5', 'A_HA_H6', 'A_HA_H7', 'A_HA_H8',
    'A_HA_H9', 'A_NA_N1', 'A_NA_N2', 'A_NA_N3', 'A_NA_N4',
    'A_NA_N5', 'A_NA_N6', 'A_NA_N7', 'A_NA_N8', 'A_NA_N9',
    'B_HA', 'B_MP', 'B_NA', 'B_NP', 'B_NS', 'B_PA', 'B_PB1', 'B_PB2'
]

def load_fasta(fasta_file):
    """Load FASTA file into dictionary."""
    fasta_dict = {}
    try:
        with open(fasta_file) as f:
            header = None
            seq = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if header:
                        fasta_dict[header] = ''.join(seq)
                    header = line[1:]
                    seq = []
                else:
                    seq.append(line)
            if header:
                fasta_dict[header] = ''.join(seq)
    except FileNotFoundError:
        print(f"Warning: {fasta_file} not found")
    return fasta_dict

def barcode_label(pathname):
    """Extract barcode from pathname."""
    m = re.search(r'(barcode\d+)', pathname)
    return m.group(1) if m else pathname

def generate_typing_table():
    """Extract HA/NA influenza subtyping from BAM files."""
    print("\n===== Generating influenza typing table =====")
    results = []
    bam_files = list(Path(REMAP_DIR).glob("*.bam"))
    
    if not bam_files:
        print(f"Warning: No BAM files found in {REMAP_DIR}")
        return
    
    for bam in sorted(bam_files):
        barcode_match = re.search(r'(barcode\d+)', bam.name)
        if not barcode_match:
            continue
        
        barcode = barcode_match.group(1)
        result = subprocess.run(['samtools', 'idxstats', str(bam)],
                              capture_output=True, text=True)
        
        best_ha, best_ha_reads = None, 0
        best_na, best_na_reads = None, 0
        influenza_b = False
        
        for line in result.stdout.strip().split('\n'):
            fields = line.split('\t')
            if len(fields) < 3:
                continue
            contig = fields[0]
            mapped = int(fields[2])
            
            if contig in ["B_HA", "B_NA"] and mapped > 0:
                influenza_b = True
            if contig.startswith("A_HA_") and mapped > best_ha_reads:
                best_ha_reads = mapped
                best_ha = contig.replace("A_HA_", "")
            if contig.startswith("A_NA_") and mapped > best_na_reads:
                best_na_reads = mapped
                best_na = contig.replace("A_NA_", "")
        
        if influenza_b:
            flu_type = "B"
        elif best_ha and best_na:
            flu_type = f"{best_ha}{best_na}"
        elif best_ha:
            flu_type = best_ha
        elif best_na:
            flu_type = best_na
        else:
            flu_type = "Unclassified"
        
        results.append({"Barcode": barcode, "Type": flu_type})
    
    if results:
        df = pd.DataFrame(results).sort_values("Barcode")
        output_file = "influenza_typing_table.csv"
        df.to_csv(output_file, index=False)
        print(f"✓ Typing table: {output_file} ({len(df)} samples)")
    else:
        print("✗ No valid barcodes found for typing")

def plot_coverage_by_barcode():
    """Plot coverage per segment by barcode (groups of 16)."""
    print("\n===== Plotting coverage by barcode =====")
    bam_files = list(Path(REMAP_DIR).glob("*.bam"))
    
    if not bam_files:
        print(f"Warning: No BAM files found in {REMAP_DIR}")
        return
    
    coverage_by_barcode = {}
    
    for bam in sorted(bam_files):
        match = re.search(r'(barcode\d+)', bam.name)
        if not match:
            continue
        barcode = match.group(1)
        coverage_by_barcode[barcode] = {}
        
        for segment in SEGMENT_LIST:
            cmd = ['samtools', 'depth', '-aa', '-r', segment, str(bam)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                coverage_by_barcode[barcode][segment] = 0
            else:
                depths = [int(line.split()[2]) for line in result.stdout.splitlines() if line]
                coverage_by_barcode[barcode][segment] = np.mean(depths) if depths else 0
    
    df = pd.DataFrame.from_dict(coverage_by_barcode, orient='index')
    df = df.reindex(columns=SEGMENT_LIST).fillna(0)
    df_filtered = df.loc[:, (df.sum() > 0)]
    df_filtered = df_filtered.sort_index()
    
    # Bar plots in chunks of 16 barcodes
    barcodes_list = df_filtered.index.tolist()
    n_barcodes = len(barcodes_list)
    n_plots = (n_barcodes + 15) // 16
    
    for plot_idx in range(n_plots):
        start_idx = plot_idx * 16
        end_idx = min((plot_idx + 1) * 16, n_barcodes)
        
        df_subset = df_filtered.iloc[start_idx:end_idx]
        
        fig, ax = plt.subplots(figsize=(14, 6))
        df_subset.plot(kind='bar', ax=ax, width=0.8)
        ax.set_title(f'Coverage by segment (barcodes {start_idx+1}-{end_idx})')
        ax.set_xlabel('Barcode')
        ax.set_ylabel('Mean depth')
        ax.legend(title='Segment', bbox_to_anchor=(1.0, 1.0), loc='upper left', fontsize=8)
        plt.tight_layout()
        plt.savefig(f"coverage_by_barcode_{plot_idx+1}.png", dpi=150, bbox_inches='tight')
        plt.close()
        print(f"✓ Coverage plot {plot_idx+1}/{n_plots}")

def plot_quality_distributions():
    """Plot MAPQ and base quality distributions."""
    print("\n===== Plotting quality distributions =====")
    bam_files = list(Path(REMAP_DIR).glob("*.bam"))
    
    if not bam_files:
        print(f"Warning: No BAM files found in {REMAP_DIR}")
        return
    
    mapq_scores, base_qualities = [], []
    max_reads_per_bam = 1000
    
    for bam in sorted(bam_files):
        proc = subprocess.Popen(['samtools', 'view', '-F', '0x4', str(bam)],
                               stdout=subprocess.PIPE, text=True)
        for i, line in enumerate(proc.stdout):
            if i >= max_reads_per_bam:
                break
            fields = line.split('\t')
            if len(fields) < 11:
                continue
            mapq_scores.append(int(fields[4]))
            qual = fields[10].rstrip()
            base_qualities.extend([ord(ch) - 33 for ch in qual])
        proc.stdout.close()
        proc.wait()
    
    if not mapq_scores or not base_qualities:
        print("✗ No quality data collected")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    
    axes[0].hist(mapq_scores, bins=range(0, 61), color="tab:blue", edgecolor="black")
    axes[0].set_title("MAPQ Distribution")
    axes[0].set_xlabel("MAPQ")
    axes[0].set_ylabel("Read Count")
    
    axes[1].hist(base_qualities, bins=np.arange(0, 51) - 0.5, color="tab:orange", edgecolor="black")
    axes[1].set_title("Base Quality Distribution")
    axes[1].set_xlabel("Phred Score")
    axes[1].set_ylabel("Base Count")
    
    plt.tight_layout()
    plt.savefig("quality_distributions.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Quality plot: quality_distributions.png")
    print(f"  MAPQ median: {np.median(mapq_scores):.1f}, mean: {np.mean(mapq_scores):.1f}")
    print(f"  Base quality median: {np.median(base_qualities):.1f}, mean: {np.mean(base_qualities):.1f}")

def plot_read_length_distribution():
    """Plot read length distribution."""
    print("\n===== Plotting read length distribution =====")
    bam_files = list(Path(REMAP_DIR).glob("*.bam"))
    
    if not bam_files:
        print(f"Warning: No BAM files found in {REMAP_DIR}")
        return
    
    read_lengths = []
    max_reads_per_bam = 15000
    
    for bam in sorted(bam_files):
        proc = subprocess.Popen(['samtools', 'view', '-F', '0x4', str(bam)],
                               stdout=subprocess.PIPE, text=True)
        for i, line in enumerate(proc.stdout):
            if i >= max_reads_per_bam:
                break
            fields = line.split('\t')
            if len(fields) < 11:
                continue
            read_lengths.append(len(fields[9]))
        proc.stdout.close()
        proc.wait()
    
    if not read_lengths:
        print("✗ No read length data collected")
        return
    
    fig, ax = plt.subplots(figsize=(10, 5))
    bins = np.linspace(0, min(max(read_lengths), 3500), 120)
    ax.hist(read_lengths, bins=bins, color="tab:blue", edgecolor="black")
    ax.set_title("Read Length Distribution (mapped reads)")
    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel("Read count")
    ax.axvline(np.median(read_lengths), color="red", linestyle="--", 
               label=f"median={np.median(read_lengths):.0f} bp")
    ax.legend()
    plt.tight_layout()
    plt.savefig("read_length_distribution.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Read length plot: read_length_distribution.png")
    print(f"  Median: {np.median(read_lengths):.0f} bp")
    print(f"  75th percentile: {np.percentile(read_lengths, 75):.0f} bp")

def plot_variant_counts():
    """Plot variant counts per barcode by segment."""
    print("\n===== Plotting variant counts =====")
    vcf_files = sorted(Path(REMAP_DIR).glob("*.filtered.vcf.gz"))
    
    if not vcf_files:
        print(f"Warning: No filtered VCF files found in {REMAP_DIR}")
        return
    
    # Simplified segment order for variant aggregation
    segment_order = ['A_HA', 'A_NA', 'A_MP', 'A_NP', 'A_NS', 'A_PA', 'A_PB1', 'A_PB2', 'B_HA', 'B_NA']
    
    def aggregate_counts(counts):
        aggregated = {}
        for seg, count in counts.items():
            if seg.startswith('A_HA_'):
                aggregated['A_HA'] = aggregated.get('A_HA', 0) + count
            elif seg.startswith('A_NA_'):
                aggregated['A_NA'] = aggregated.get('A_NA', 0) + count
            elif seg.startswith('B_HA'):
                aggregated['B_HA'] = aggregated.get('B_HA', 0) + count
            elif seg.startswith('B_NA'):
                aggregated['B_NA'] = aggregated.get('B_NA', 0) + count
            else:
                aggregated[seg] = count
        return aggregated
    
    variant_counts = {}
    for vcf in vcf_files:
        sample = barcode_label(vcf.name)
        
        result = subprocess.run(['bcftools', 'query', '-f', '%CHROM\n', str(vcf)],
                              capture_output=True, text=True)
        if result.returncode != 0:
            continue
        
        counts = Counter(result.stdout.splitlines())
        aggregated = aggregate_counts(counts)
        variant_counts[sample] = {seg: aggregated.get(seg, 0) for seg in segment_order}
    
    if not variant_counts:
        print("✗ No variant data collected")
        return
    
    variant_df = pd.DataFrame.from_dict(variant_counts, orient='index')
    variant_df = variant_df.reindex(columns=segment_order).fillna(0).sort_index()
    
    # Save as CSV
    variant_df.to_csv("variant_counts.csv")
    
    # Plot in groups of 16
    n_samples = len(variant_df)
    n_plots = (n_samples + 15) // 16
    for plot_idx in range(n_plots):
        start_idx = plot_idx * 16
        end_idx = min((plot_idx + 1) * 16, n_samples)
        subset = variant_df.iloc[start_idx:end_idx]
        
        fig, ax = plt.subplots(figsize=(14, 6))
        subset.plot(kind='bar', stacked=True, ax=ax, width=0.8, colormap=mpl.cm.get_cmap('tab20'))
        ax.set_title(f'Variant counts by segment (barcodes {start_idx+1}-{end_idx})')
        ax.set_xlabel('Barcode')
        ax.set_ylabel('Variant count')
        ax.legend(title='Segment', bbox_to_anchor=(1.0, 1.0), loc='upper left', fontsize=9)
        plt.tight_layout()
        plt.savefig(f"variant_counts_{plot_idx+1}.png", dpi=150, bbox_inches='tight')
        plt.close()
        print(f"✓ Variant plot {plot_idx+1}/{n_plots}")

def plot_consensus_identity():
    """Plot consensus sequence identity to reference."""
    print("\n===== Plotting consensus identity =====")
    
    if not os.path.exists(REFERENCE_FASTA):
        print(f"Warning: Reference FASTA {REFERENCE_FASTA} not found")
        return
    
    ref_seqs = load_fasta(REFERENCE_FASTA)
    consensus_files = sorted(Path(REMAP_DIR).glob("*.consensus.fasta"))
    
    if not consensus_files:
        print(f"Warning: No consensus files found in {REMAP_DIR}")
        return
    
    identity_by_sample = {}
    for fasta in consensus_files:
        sample = barcode_label(fasta.name)
        seqs = load_fasta(str(fasta))
        
        total_bases = 0
        total_matches = 0
        for contig, consensus_seq in seqs.items():
            if contig not in ref_seqs:
                continue
            ref_seq = ref_seqs[contig]
            pair_len = min(len(ref_seq), len(consensus_seq))
            if pair_len == 0:
                continue
            total_bases += pair_len
            total_matches += sum(1 for a, b in zip(ref_seq[:pair_len], consensus_seq[:pair_len]) if a == b)
        
        if total_bases:
            identity_by_sample[sample] = total_matches / total_bases * 100
    
    if not identity_by_sample:
        print("✗ No consensus identity data collected")
        return
    
    identity_series = pd.Series(identity_by_sample).sort_index()
    
    fig, ax = plt.subplots(figsize=(14, 5))
    identity_series.plot(kind='bar', ax=ax, color='tab:blue')
    ax.set_title('Consensus Sequence Identity to Reference')
    ax.set_xlabel('Barcode')
    ax.set_ylabel('Identity (%)')
    ax.set_ylim(max(0, identity_series.min() - 2), 100.5)
    plt.tight_layout()
    plt.savefig("consensus_identity.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Consensus identity plot: consensus_identity.png")
    print(f"  Mean identity: {identity_series.mean():.2f}%")

def plot_coverage_heatmap():
    """Plot coverage heatmap for segments."""
    print("\n===== Plotting coverage heatmap =====")
    bam_files = list(Path(REMAP_DIR).glob("*.bam"))
    
    if not bam_files:
        print(f"Warning: No BAM files found in {REMAP_DIR}")
        return
    
    coverage_by_barcode = {}
    for bam in sorted(bam_files):
        match = re.search(r'(barcode\d+)', bam.name)
        if not match:
            continue
        barcode = match.group(1)
        coverage_by_barcode[barcode] = {}
        
        for segment in SEGMENT_LIST:
            cmd = ['samtools', 'depth', '-aa', '-r', segment, str(bam)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                coverage_by_barcode[barcode][segment] = 0
            else:
                depths = [int(line.split()[2]) for line in result.stdout.splitlines() if line]
                coverage_by_barcode[barcode][segment] = np.mean(depths) if depths else 0
    
    df = pd.DataFrame.from_dict(coverage_by_barcode, orient='index')
    df = df.reindex(columns=SEGMENT_LIST).fillna(0)
    df_filtered = df.loc[:, (df.sum() > 0)]
    
    # Heatmaps in chunks of 16 segments
    segments_list = df_filtered.columns.tolist()
    n_segments = len(segments_list)
    n_plots = (n_segments + 15) // 16
    
    for plot_idx in range(n_plots):
        start_idx = plot_idx * 16
        end_idx = min((plot_idx + 1) * 16, n_segments)
        chunk = segments_list[start_idx:end_idx]
        
        heat_data = np.log10(df_filtered[chunk] + 1).T
        
        fig, ax = plt.subplots(figsize=(14, 8))
        im = ax.imshow(heat_data, aspect='auto', cmap='magma', interpolation='nearest')
        
        ax.set_xticks(np.arange(len(heat_data.columns)))
        ax.set_xticklabels(heat_data.columns, rotation=90, fontsize=8)
        ax.set_yticks(np.arange(len(chunk)))
        ax.set_yticklabels(chunk, fontsize=9)
        
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('log10(mean depth + 1)')
        
        ax.set_title(f'Coverage heatmap (segments {start_idx+1}-{end_idx})')
        ax.set_xlabel('Barcode')
        ax.set_ylabel('Segment')
        
        plt.tight_layout()
        plt.savefig(f"coverage_heatmap_{plot_idx+1}.png", dpi=150, bbox_inches='tight')
        plt.close()
        print(f"✓ Heatmap {plot_idx+1}/{n_plots}")

def main():
    """Run all analysis and plotting."""
    print("=" * 60)
    print("INFLUENZA PIPELINE: ANALYSIS AND VISUALIZATION")
    print("=" * 60)
    
    try:
        generate_typing_table()
        plot_coverage_by_barcode()
        plot_quality_distributions()
        plot_read_length_distribution()
        plot_variant_counts()
        plot_consensus_identity()
        plot_coverage_heatmap()
        
        print("\n" + "=" * 60)
        print("✓ All analysis complete!")
        print("=" * 60)
        print("\nGenerated files:")
        print("  - influenza_typing_table.csv")
        print("  - coverage_by_barcode_*.png")
        print("  - quality_distributions.png")
        print("  - read_length_distribution.png")
        print("  - variant_counts.csv and variant_counts_*.png")
        print("  - consensus_identity.png")
        print("  - coverage_heatmap_*.png")
        
    except Exception as e:
        print(f"\n✗ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == '__main__':
    exit(main())
