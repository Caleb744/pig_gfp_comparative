#!/usr/bin/env python3
"""
per_barcode_read_distribution.py
Usage:
    conda activate /cluster/pixstor/schnabelr-drii/WORKING/cs744/CONDA_ENV/barcode_reads
    python per_barcode_read_distribution.py \
      --bam possorted_genome_bam.bam \
      --gtf genes.gtf \
      --out per_barcode_readdist.tsv \
      --threads 4

Outputs TSV with columns:
CB, total_reads, exonic_reads, intronic_reads, intergenic_reads, mito_reads
"""

import argparse
import pysam
from intervaltree import IntervalTree
import gzip
import sys
from collections import defaultdict

def build_gene_intervals(gtf_path):
    # Build two structures:
    #  - exon_trees[chr] = IntervalTree of exons (for exonic overlap)
    #  - gene_trees[chr] = IntervalTree of genes (for intronic overlap detection)
    exon_trees = defaultdict(IntervalTree)
    gene_trees = defaultdict(IntervalTree)

    openf = gzip.open if gtf_path.endswith('.gz') else open
    with openf(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            chrom, source, feat, start, end, score, strand, frame, attr = line.strip().split('\t')
            s = int(start) - 1  # convert to 0-based
            e = int(end)
            if feat == 'exon':
                exon_trees[chrom].addi(s, e, None)
            elif feat == 'gene':
                gene_trees[chrom].addi(s, e, None)

    return exon_trees, gene_trees

def classify_read(chrom, start, end, exon_trees, gene_trees):
    # Return "exonic" / "intronic" / "intergenic"
    if chrom not in exon_trees and chrom not in gene_trees:
        return "intergenic"
    if chrom in exon_trees and exon_trees[chrom].overlap(start, end):
        return "exonic"
    if chrom in gene_trees and gene_trees[chrom].overlap(start, end):
        return "intronic"
    return "intergenic"

def main(args):
    exon_trees, gene_trees = build_gene_intervals(args.gtf)
    bam = pysam.AlignmentFile(args.bam, "rb")
    counts = defaultdict(lambda: {'total':0, 'exonic':0, 'intronic':0, 'intergenic':0, 'mito':0})

    # iterate reads, skip secondary/supplementary/unmapped
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped: continue
        if read.is_secondary or read.is_supplementary: continue

        # find CB tag for CellRanger-style BAM; fallback to "CR" if needed
        try:
            cb = read.get_tag('CB')
        except KeyError:
            try:
                cb = read.get_tag('CR')  # rarely used
            except KeyError:
                continue

        counts[cb]['total'] += 1

        chrom = bam.get_reference_name(read.reference_id)
        start = read.reference_start
        end = read.reference_end

        if chrom in ('chrM', 'MT', 'M'):
            counts[cb]['mito'] += 1

        cls = classify_read(chrom, start, end, exon_trees, gene_trees)
        counts[cb][cls] += 1

    # write out
    with open(args.out, 'w') as out:
        out.write("CB\ttotal_reads\texonic_reads\tintronic_reads\tintergenic_reads\tmito_reads\n")
        for cb, d in counts.items():
            out.write(f"{cb}\t{d['total']}\t{d['exonic']}\t{d['intronic']}\t{d['intergenic']}\t{d['mito']}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', required=True)
    parser.add_argument('--gtf', required=True, help="GTF/GFF of genes (exon + gene features)")
    parser.add_argument('--out', required=True)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    main(args)
