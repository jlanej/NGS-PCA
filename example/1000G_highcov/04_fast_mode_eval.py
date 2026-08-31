#!/usr/bin/env python3
"""Evaluate the mosdepth fast-mode comparison.

Compares two NGS-PCA runs - one from mosdepth as configured, one from
mosdepth --fast-mode - and the paired mosdepth wall times recorded by
01_download_and_mosdepth.sh. Writes:

  pc_correlation.tsv    per-PC Pearson r between the runs, and the ratio of
                        singular values
  timing_summary.tsv    per-mode wall-time statistics and paired speedups
  qc_concordance.tsv    per-column agreement of the two sample_qc.tsv tables,
                        when both were produced (see the README) - this is
                        where MTDNA_CN and the coverage ratios are compared
  fast_mode_report.md   the headline numbers and their caveats
  fast_mode_summary.png correlations, scatters, runtime boxplot, speedup
                        histogram (only when matplotlib is available)
  fast_mode_qc.png      QC column agreement and the MTDNA_CN scatter (only
                        with matplotlib and both QC tables)

The sign of a singular vector is arbitrary, so correlations are reported as
computed and evaluated as |r|. Needs only the standard library; matplotlib is
optional and its absence skips the figures, not the evaluation.
"""

import argparse
import glob
import math
import os
import statistics
import sys


def read_pcs(path):
    """svd.pcs.txt -> (ordered sample list, {sample: [pc values]}, pc count)"""
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if not header or header[0] != "SAMPLE":
            sys.exit(f"ERROR: {path} does not look like svd.pcs.txt (header: {header[:3]}...)")
        n_pcs = len(header) - 1
        samples, values = [], {}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != n_pcs + 1:
                sys.exit(f"ERROR: {path} row for {fields[0]} has {len(fields) - 1} PCs, header has {n_pcs}")
            samples.append(fields[0])
            values[fields[0]] = [float(v) for v in fields[1:]]
    return samples, values, n_pcs


def read_singular_values(path):
    with open(path) as handle:
        handle.readline()
        return [float(line.split("\t")[1]) for line in handle if line.strip()]


def pearson(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def read_timing(timing_dir):
    """Timing rows -> {sample: {mode: {"wall_s": int, "first_mode": str}}}, versions seen"""
    rows, versions = {}, set()
    for path in sorted(glob.glob(os.path.join(timing_dir, "*.tsv"))):
        with open(path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue
                sample, mode, first_mode, wall_s = fields[0], fields[1], fields[2], fields[3]
                rows.setdefault(sample, {})[mode] = {
                    "wall_s": float(wall_s),
                    "first_mode": first_mode,
                }
                if len(fields) >= 5:
                    versions.add(fields[4])
    return rows, versions


def summarise(values):
    if not values:
        return {"n": 0}
    ordered = sorted(values)
    return {
        "n": len(values),
        "median": statistics.median(ordered),
        "mean": sum(ordered) / len(ordered),   # not statistics.fmean: that needs python 3.8
        "sd": statistics.pstdev(ordered) if len(ordered) > 1 else 0.0,
        "min": ordered[0],
        "max": ordered[-1],
    }


def same_file_content(a, b):
    try:
        with open(a) as fa, open(b) as fb:
            return fa.read() == fb.read()
    except OSError:
        return False


def read_qc(path):
    """sample_qc.tsv -> (column names after SAMPLE_ID, {sample: {column: text}})"""
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = {}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            rows[fields[0]] = dict(zip(header[1:], fields[1:]))
    return header[1:], rows


def qc_concordance(normal_qc, fast_qc):
    """Per-column agreement between the two QC tables.

    Numeric columns get Pearson r and the median signed relative difference
    (fast minus normal, as a percent of normal) - correlation alone would
    miss a uniform bias, which is exactly what skipping mate-overlap
    correction produces. Non-numeric columns are skipped, except INFERRED_SEX,
    which is reported as an agreement count.
    """
    normal_columns, normal_rows = normal_qc
    fast_columns, fast_rows = fast_qc
    shared = [s for s in normal_rows if s in fast_rows]
    numeric, sex = [], None
    for column in [c for c in normal_columns if c in set(fast_columns)]:
        pairs = []
        for sample in shared:
            try:
                pairs.append((float(normal_rows[sample][column]),
                              float(fast_rows[sample][column])))
            except (KeyError, ValueError):
                pass
        if len(pairs) >= 3 and len(pairs) >= 0.9 * len(shared):
            xs = [p[0] for p in pairs]
            ys = [p[1] for p in pairs]
            rel = [100.0 * (y - x) / abs(x) for x, y in pairs if x != 0]
            numeric.append({"column": column, "n": len(pairs), "r": pearson(xs, ys),
                            "median_rel_diff_pct": statistics.median(rel) if rel else float("nan")})
        elif column == "INFERRED_SEX":
            agree = sum(1 for s in shared
                        if normal_rows[s].get(column) == fast_rows[s].get(column))
            sex = (agree, len(shared))
    return len(shared), numeric, sex


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--normal", required=True, help="NGS-PCA output dir from standard mosdepth")
    parser.add_argument("--fast", required=True, help="NGS-PCA output dir from mosdepth --fast-mode")
    parser.add_argument("--timing", default=None, help="directory of per-sample timing rows")
    parser.add_argument("--qc-normal", default=None, help="sample_qc.tsv from the standard tree")
    parser.add_argument("--qc-fast", default=None, help="sample_qc.tsv from the fast tree")
    parser.add_argument("--out", required=True, help="directory to write the evaluation into")
    parser.add_argument("--pcs", type=int, default=0,
                        help="leading PCs to evaluate (default 0 = every PC both runs computed)")
    args = parser.parse_args()

    for side, directory in (("normal", args.normal), ("fast", args.fast)):
        if not os.path.isfile(os.path.join(directory, "svd.pcs.txt")):
            sys.exit(f"ERROR: {directory} has no svd.pcs.txt - has the {side} NGS-PCA run finished?")
    os.makedirs(args.out, exist_ok=True)

    normal_samples, normal_pcs, normal_n = read_pcs(os.path.join(args.normal, "svd.pcs.txt"))
    fast_samples, fast_pcs, fast_n = read_pcs(os.path.join(args.fast, "svd.pcs.txt"))
    warnings = []

    shared = [s for s in normal_samples if s in fast_pcs]
    if len(shared) < 2:
        sys.exit(f"ERROR: only {len(shared)} samples are shared between the two runs - nothing to correlate")
    if len(shared) != len(normal_samples) or len(shared) != len(fast_samples):
        warnings.append(f"sample sets differ: {len(normal_samples)} normal, {len(fast_samples)} fast, "
                        f"{len(shared)} shared - correlations use the shared samples only")

    # The two runs select bins from their own first file; the selection depends
    # only on bin positions, which fast mode does not change, so anything but
    # identical here means the trees were not built from the same cohort.
    bins_identical = same_file_content(os.path.join(args.normal, "svd.bins.txt"),
                                       os.path.join(args.fast, "svd.bins.txt"))
    if not bins_identical:
        warnings.append("svd.bins.txt differs between the runs - the two PCAs did not use the same "
                        "bins, so PC correlations conflate fast mode with bin selection")

    k = min(args.pcs, normal_n, fast_n) if args.pcs > 0 else min(normal_n, fast_n)
    sv_normal = read_singular_values(os.path.join(args.normal, "svd.singularvalues.txt"))
    sv_fast = read_singular_values(os.path.join(args.fast, "svd.singularvalues.txt"))

    correlations = []
    for pc in range(k):
        xs = [normal_pcs[s][pc] for s in shared]
        ys = [fast_pcs[s][pc] for s in shared]
        r = pearson(xs, ys)
        ratio = (sv_fast[pc] / sv_normal[pc]) if pc < min(len(sv_normal), len(sv_fast)) and sv_normal[pc] else float("nan")
        correlations.append({"pc": pc + 1, "r": r, "abs_r": abs(r), "sv_ratio": ratio})

    with open(os.path.join(args.out, "pc_correlation.tsv"), "w") as out:
        out.write("PC\tpearson_r\tabs_r\tsingular_value_ratio_fast_over_normal\n")
        for row in correlations:
            out.write(f"{row['pc']}\t{row['r']:.6f}\t{row['abs_r']:.6f}\t{row['sv_ratio']:.6f}\n")

    finite = [c["abs_r"] for c in correlations if not math.isnan(c["abs_r"])]
    min_abs_r = min(finite) if finite else float("nan")
    weakest = min((c for c in correlations if not math.isnan(c["abs_r"])),
                  key=lambda c: c["abs_r"], default=None)

    # ── Timing ────────────────────────────────────────────────────────────────
    timing, versions = ({}, set())
    if args.timing and os.path.isdir(args.timing):
        timing, versions = read_timing(args.timing)
    elif args.timing:
        warnings.append(f"timing directory not found: {args.timing} - runtime comparison skipped")

    normal_walls = [t["normal"]["wall_s"] for t in timing.values() if "normal" in t]
    fast_walls = [t["fast"]["wall_s"] for t in timing.values() if "fast" in t]
    paired = {s: t for s, t in timing.items() if "normal" in t and "fast" in t and t["fast"]["wall_s"] > 0}
    speedups = [t["normal"]["wall_s"] / t["fast"]["wall_s"] for t in paired.values()]
    by_order = {"normal": [], "fast": []}
    for t in paired.values():
        by_order[t["normal"]["first_mode"]].append(t["normal"]["wall_s"] / t["fast"]["wall_s"])

    # ── QC concordance ────────────────────────────────────────────────────────
    qc_shared, qc_numeric, qc_sex = 0, [], None
    qc_note = None
    if args.qc_normal and args.qc_fast:
        if os.path.isfile(args.qc_normal) and os.path.isfile(args.qc_fast):
            qc_shared, qc_numeric, qc_sex = qc_concordance(read_qc(args.qc_normal),
                                                           read_qc(args.qc_fast))
            with open(os.path.join(args.out, "qc_concordance.tsv"), "w") as out:
                out.write("column\tn\tpearson_r\tmedian_rel_diff_pct\n")
                for row in qc_numeric:
                    out.write(f"{row['column']}\t{row['n']}\t{row['r']:.6f}\t"
                              f"{row['median_rel_diff_pct']:.2f}\n")
        else:
            qc_note = ("QC tables not found - run 03a and 03 for both trees (see the README) to "
                       "compare depth-derived phenotypes such as MTDNA_CN")

    normal_stats, fast_stats, speedup_stats = summarise(normal_walls), summarise(fast_walls), summarise(speedups)
    with open(os.path.join(args.out, "timing_summary.tsv"), "w") as out:
        out.write("metric\tn\tmedian\tmean\tsd\tmin\tmax\n")
        for name, stats in (("normal_wall_s", normal_stats), ("fast_wall_s", fast_stats),
                            ("paired_speedup", speedup_stats)):
            if stats["n"]:
                out.write(f"{name}\t{stats['n']}\t{stats['median']:.2f}\t{stats['mean']:.2f}\t"
                          f"{stats['sd']:.2f}\t{stats['min']:.2f}\t{stats['max']:.2f}\n")
            else:
                out.write(f"{name}\t0\t\t\t\t\t\n")

    # ── Report ────────────────────────────────────────────────────────────────
    lines = ["# mosdepth fast-mode evaluation", ""]
    lines.append(f"{len(shared)} shared samples; first {k} PCs evaluated; "
                 f"mosdepth version(s) in timing rows: {', '.join(sorted(versions)) or 'not recorded'}.")
    lines.append("")
    lines.append(f"- Minimum |r| over PC1-PC{k}: **{min_abs_r:.6f}**"
                 + (f" (weakest: PC{weakest['pc']})" if weakest else ""))
    if speedup_stats["n"]:
        lines.append(f"- Paired samples with both timings: {speedup_stats['n']}")
        lines.append(f"- Median wall time: normal {normal_stats['median']:.0f} s, "
                     f"fast {fast_stats['median']:.0f} s")
        lines.append(f"- Median per-sample speedup (normal/fast): **{speedup_stats['median']:.2f}x** "
                     f"(range {speedup_stats['min']:.2f}-{speedup_stats['max']:.2f})")
        if by_order["normal"] and by_order["fast"]:
            lines.append(f"- Order check - median speedup when normal ran first: "
                         f"{statistics.median(by_order['normal']):.2f}x; when fast ran first: "
                         f"{statistics.median(by_order['fast']):.2f}x. A large gap would mean cache "
                         f"warming, not fast mode, is being measured.")
    else:
        lines.append("- No paired timing rows found; runtime comparison skipped.")
    if qc_numeric:
        lines.append("")
        lines.append(f"## QC concordance ({qc_shared} shared samples)")
        lines.append("")
        mtdna = next((row for row in qc_numeric if row["column"] == "MTDNA_CN"), None)
        if mtdna:
            lines.append(f"- MTDNA_CN: r = **{mtdna['r']:.6f}**, median relative difference "
                         f"fast vs normal = **{mtdna['median_rel_diff_pct']:+.2f}%**. A nonzero "
                         f"shift here is expected in kind if not in size: skipping mate-overlap "
                         f"correction inflates chrM and the nuclear genome by different amounts "
                         f"when their insert sizes differ. Correlation, not the shift, decides "
                         f"whether fast-mode MTDNA_CN ranks samples the same way.")
        if qc_sex:
            lines.append(f"- INFERRED_SEX agreement: {qc_sex[0]}/{qc_sex[1]}")
        lines.append(f"- Every numeric column is in qc_concordance.tsv ({len(qc_numeric)} compared).")
    elif qc_note:
        lines.append("")
        lines.append(f"- {qc_note}.")
    lines.append("")
    lines.append("Fast mode skips CIGAR-aware depth and mate-overlap correction. Both act mostly "
                 "as per-sample multiplicative shifts, which the log2 fold change against each "
                 "sample's own median absorbs - so near-1 correlations are the expectation, and "
                 "a PC that breaks from it is worth inspecting alongside its singular-value ratio: "
                 "adjacent PCs with near-equal singular values can rotate into each other without "
                 "any real disagreement between the runs.")
    for warning in warnings:
        lines.append("")
        lines.append(f"**WARNING:** {warning}")
    with open(os.path.join(args.out, "fast_mode_report.md"), "w") as out:
        out.write("\n".join(lines) + "\n")

    # ── Figure (optional) ─────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available - skipping fast_mode_summary.png")
        plt = None

    if plt is not None:
        fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
        pcs_axis = [c["pc"] for c in correlations]

        axes[0][0].bar(pcs_axis, [c["abs_r"] for c in correlations], color="#33689E")
        axes[0][0].set_ylim(min(0.95, min_abs_r * 0.999) if finite else 0, 1.001)
        axes[0][0].set_xlabel("PC")
        axes[0][0].set_ylabel("|Pearson r|")
        axes[0][0].set_title(f"PC correlation, normal vs fast (min {min_abs_r:.4f})")

        for panel, pc in ((axes[0][1], 0), (axes[0][2], 1)):
            if pc >= k:
                panel.axis("off")
                continue
            xs = [normal_pcs[s][pc] for s in shared]
            sign = -1.0 if correlations[pc]["r"] < 0 else 1.0
            ys = [sign * fast_pcs[s][pc] for s in shared]
            panel.plot(xs, ys, ".", markersize=3, alpha=0.6, color="#33689E")
            panel.set_xlabel(f"PC{pc + 1} normal")
            panel.set_ylabel(f"PC{pc + 1} fast (sign-aligned)")
            panel.set_title(f"PC{pc + 1}: r = {correlations[pc]['r']:.5f}")

        if normal_walls and fast_walls:
            axes[1][0].boxplot([normal_walls, fast_walls])
            axes[1][0].set_xticks([1, 2])
            axes[1][0].set_xticklabels(["normal", "fast"])
            axes[1][0].set_ylabel("mosdepth wall time (s)")
            axes[1][0].set_title(f"Runtime across {len(paired)} paired samples")
        else:
            axes[1][0].axis("off")

        if speedups:
            axes[1][1].hist(speedups, bins=30, color="#33689E")
            axes[1][1].set_xlabel("per-sample speedup (normal / fast)")
            axes[1][1].set_ylabel("samples")
            axes[1][1].set_title(f"Median speedup {speedup_stats['median']:.2f}x")
        else:
            axes[1][1].axis("off")

        ratios = [c["sv_ratio"] for c in correlations]
        axes[1][2].plot(pcs_axis, ratios, "o-", markersize=3, color="#33689E")
        axes[1][2].axhline(1.0, color="grey", linewidth=0.8)
        axes[1][2].set_xlabel("PC")
        axes[1][2].set_ylabel("singular value ratio (fast / normal)")
        axes[1][2].set_title("Scale agreement per PC")

        fig.suptitle("mosdepth fast-mode comparison", fontsize=14)
        fig.savefig(os.path.join(args.out, "fast_mode_summary.png"), dpi=150)
        print(f"wrote {os.path.join(args.out, 'fast_mode_summary.png')}")

    if plt is not None and qc_numeric:
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
        names = [row["column"] for row in qc_numeric]

        axes[0].bar(range(len(qc_numeric)), [abs(row["r"]) for row in qc_numeric], color="#33689E")
        axes[0].set_xticks(range(len(qc_numeric)))
        axes[0].set_xticklabels(names, rotation=60, ha="right", fontsize=7)
        axes[0].set_ylabel("|Pearson r|")
        axes[0].set_title("QC column correlation, normal vs fast")

        scatter = next((row for row in qc_numeric if row["column"] == "MTDNA_CN"), qc_numeric[0])
        normal_columns, normal_rows = read_qc(args.qc_normal)
        fast_columns, fast_rows = read_qc(args.qc_fast)
        pairs = []
        for sample in normal_rows:
            if sample in fast_rows:
                try:
                    pairs.append((float(normal_rows[sample][scatter["column"]]),
                                  float(fast_rows[sample][scatter["column"]])))
                except (KeyError, ValueError):
                    pass
        xs = [p[0] for p in pairs]
        ys = [p[1] for p in pairs]
        axes[1].plot(xs, ys, ".", markersize=4, alpha=0.6, color="#33689E")
        lo, hi = min(xs + ys), max(xs + ys)
        axes[1].plot([lo, hi], [lo, hi], color="grey", linewidth=0.8)
        axes[1].set_xlabel(f"{scatter['column']} normal")
        axes[1].set_ylabel(f"{scatter['column']} fast")
        axes[1].set_title(f"{scatter['column']}: r = {scatter['r']:.5f}, "
                          f"median shift {scatter['median_rel_diff_pct']:+.2f}%")

        axes[2].bar(range(len(qc_numeric)),
                    [row["median_rel_diff_pct"] for row in qc_numeric], color="#33689E")
        axes[2].axhline(0.0, color="grey", linewidth=0.8)
        axes[2].set_xticks(range(len(qc_numeric)))
        axes[2].set_xticklabels(names, rotation=60, ha="right", fontsize=7)
        axes[2].set_ylabel("median relative difference (%)")
        axes[2].set_title("Systematic shift, fast vs normal")

        fig.suptitle("sample_qc.tsv concordance", fontsize=13)
        fig.savefig(os.path.join(args.out, "fast_mode_qc.png"), dpi=150)
        print(f"wrote {os.path.join(args.out, 'fast_mode_qc.png')}")

    print(f"wrote {os.path.join(args.out, 'pc_correlation.tsv')}")
    print(f"wrote {os.path.join(args.out, 'timing_summary.tsv')}")
    if qc_numeric:
        print(f"wrote {os.path.join(args.out, 'qc_concordance.tsv')}")
    print(f"wrote {os.path.join(args.out, 'fast_mode_report.md')}")
    print(f"min |r| over PC1-PC{k}: {min_abs_r:.6f}"
          + (f"; median speedup {speedup_stats['median']:.2f}x" if speedup_stats["n"] else ""))


if __name__ == "__main__":
    main()
