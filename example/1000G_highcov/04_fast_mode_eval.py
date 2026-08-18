#!/usr/bin/env python3
"""Evaluate the mosdepth fast-mode comparison.

Compares two NGS-PCA runs - one from mosdepth as configured, one from
mosdepth --fast-mode - and the paired mosdepth wall times recorded by
01_download_and_mosdepth.sh. Writes:

  pc_correlation.tsv    per-PC Pearson r between the runs, and the ratio of
                        singular values
  timing_summary.tsv    per-mode wall-time statistics and paired speedups
  fast_mode_report.md   the headline numbers and their caveats
  fast_mode_summary.png correlations, scatters, runtime boxplot, speedup
                        histogram (only when matplotlib is available)

The sign of a singular vector is arbitrary, so correlations are reported as
computed and evaluated as |r|. Needs only the standard library; matplotlib is
optional and its absence skips the figure, not the evaluation.
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
        "mean": statistics.fmean(ordered),
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


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--normal", required=True, help="NGS-PCA output dir from standard mosdepth")
    parser.add_argument("--fast", required=True, help="NGS-PCA output dir from mosdepth --fast-mode")
    parser.add_argument("--timing", default=None, help="directory of per-sample timing rows")
    parser.add_argument("--out", required=True, help="directory to write the evaluation into")
    parser.add_argument("--pcs", type=int, default=20, help="leading PCs to evaluate (default 20)")
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

    k = min(args.pcs, normal_n, fast_n)
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

    print(f"wrote {os.path.join(args.out, 'pc_correlation.tsv')}")
    print(f"wrote {os.path.join(args.out, 'timing_summary.tsv')}")
    print(f"wrote {os.path.join(args.out, 'fast_mode_report.md')}")
    print(f"min |r| over PC1-PC{k}: {min_abs_r:.6f}"
          + (f"; median speedup {speedup_stats['median']:.2f}x" if speedup_stats["n"] else ""))


if __name__ == "__main__":
    main()
