#!/usr/bin/env python3
"""
Compare UNIFORM vs GAUSSIAN random matrix distributions for randomized SVD PCA.

Reads svd.pcs.txt and svd.singularvalues.txt from two output directories
(one run with --distribution UNIFORM, one with --distribution GAUSSIAN) and
writes a TSV summary with per-PC descriptive statistics and Pearson/Spearman
correlations. The absolute correlation is reported because PCs can differ by
a sign flip (which is mathematically equivalent).

Usage:
    python compare_distributions.py <uniform_dir> <gaussian_dir> <output_file>
"""

import sys
import os
import math

def read_pcs(path):
    """Return {sample: [pc1, pc2, ...]} and ordered list of pc names."""
    samples = {}
    pc_names = []
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        pc_names = header[1:]
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            samples[parts[0]] = [float(x) for x in parts[1:]]
    return samples, pc_names

def read_singular_values(path):
    """Return {pc_number: singular_value}."""
    sv = {}
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            sv[int(parts[0])] = float(parts[1])
    return sv

def mean(vals):
    return sum(vals) / len(vals)

def median(vals):
    s = sorted(vals)
    n = len(s)
    mid = n // 2
    return (s[mid - 1] + s[mid]) / 2.0 if n % 2 == 0 else s[mid]

def sd(vals):
    m = mean(vals)
    variance = sum((x - m) ** 2 for x in vals) / (len(vals) - 1)
    return math.sqrt(variance)

def pearson(x, y):
    n = len(x)
    mx, my = mean(x), mean(y)
    num = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    den = math.sqrt(sum((xi - mx) ** 2 for xi in x) * sum((yi - my) ** 2 for yi in y))
    return num / den if den != 0 else float('nan')

def rank(vals):
    """Return rank list (1-based, average for ties)."""
    indexed = sorted(enumerate(vals), key=lambda t: t[1])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(indexed):
        j = i
        while j + 1 < len(indexed) and indexed[j + 1][1] == indexed[i][1]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks[indexed[k][0]] = avg_rank
        i = j + 1
    return ranks

def spearman(x, y):
    return pearson(rank(x), rank(y))

def main():
    if len(sys.argv) != 4:
        print("Usage: compare_distributions.py <uniform_dir> <gaussian_dir> <output_file>",
              file=sys.stderr)
        sys.exit(1)

    uniform_dir = sys.argv[1]
    gaussian_dir = sys.argv[2]
    output_file = sys.argv[3]

    uniform_pcs, pc_names = read_pcs(os.path.join(uniform_dir, "svd.pcs.txt"))
    gaussian_pcs, _ = read_pcs(os.path.join(gaussian_dir, "svd.pcs.txt"))

    uniform_sv = read_singular_values(os.path.join(uniform_dir, "svd.singularvalues.txt"))
    gaussian_sv = read_singular_values(os.path.join(gaussian_dir, "svd.singularvalues.txt"))

    # Ensure same sample order
    samples = sorted(uniform_pcs.keys())
    if sorted(gaussian_pcs.keys()) != samples:
        raise ValueError("Sample sets differ between the two runs")

    lines = []
    lines.append("# UNIFORM vs GAUSSIAN random matrix distribution comparison")
    lines.append(f"# n_samples={len(samples)}")
    lines.append("# Correlations use absolute value to account for arbitrary PC sign flips")
    lines.append("#")
    header_cols = [
        "PC",
        "UNIFORM_SINGULAR_VALUE", "GAUSSIAN_SINGULAR_VALUE",
        "UNIFORM_MEAN", "UNIFORM_MEDIAN", "UNIFORM_SD",
        "GAUSSIAN_MEAN", "GAUSSIAN_MEDIAN", "GAUSSIAN_SD",
        "ABS_PEARSON_R", "ABS_SPEARMAN_R",
    ]
    lines.append("\t".join(header_cols))

    for pc_idx, pc_name in enumerate(pc_names):
        u_vals = [uniform_pcs[s][pc_idx] for s in samples]
        g_vals = [gaussian_pcs[s][pc_idx] for s in samples]
        pc_num = pc_idx + 1

        pr = abs(pearson(u_vals, g_vals))
        sr = abs(spearman(u_vals, g_vals))

        u_sv = uniform_sv.get(pc_num, float('nan'))
        g_sv = gaussian_sv.get(pc_num, float('nan'))

        row = [
            pc_name,
            f"{u_sv:.6f}", f"{g_sv:.6f}",
            f"{mean(u_vals):.6f}", f"{median(u_vals):.6f}", f"{sd(u_vals):.6f}",
            f"{mean(g_vals):.6f}", f"{median(g_vals):.6f}", f"{sd(g_vals):.6f}",
            f"{pr:.6f}", f"{sr:.6f}",
        ]
        lines.append("\t".join(row))

    with open(output_file, 'w') as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"Written: {output_file}")
    # Print a human-readable summary to stdout
    for line in lines:
        print(line)

if __name__ == "__main__":
    main()
