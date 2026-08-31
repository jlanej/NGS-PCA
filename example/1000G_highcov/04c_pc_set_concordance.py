#!/usr/bin/env python3
"""Concordance of two runs' PC sets as subspaces, not as ordered lists.

Run by 04b_fast_mode_report.sh inside the eval container (numpy arrives with
matplotlib there); the report assembler embeds what this writes. The two runs
are any pair of svd.pcs.txt tables - fast vs normal, or the seed-control
calibration pair (the same tree under two -randomSeed values), relabeled via
--label-a/--label-b.

Matched-index correlation is the strictest comparison and the wrong one deep
in the spectrum: near-equal singular values let adjacent PCs rotate into each
other (or swap rank outright) without the span - the thing downstream
regression actually uses - changing at all. This measures the span.

Inputs: the two svd.pcs.txt tables (samples x PCs), joined by sample.
Outputs, written to --out-dir:

  subspace_angles.tsv     for nested leading-k subspaces: canonical
                          correlations' minimum (cos of the largest principal
                          angle), that angle in degrees, the mean cos, and
                          how many of the k cosines fall below 0.99 and 0.90
                          - one straddled boundary and a broad divergence
                          both make a large angle, and only the counts tell
                          them apart. cos ~ 1 at every k means the
                          k-dimensional leading subspaces are interchangeable.
  pc_containment.tsv      per run-A PC: R^2 of its projection onto the span
                          of ALL run-B PCs (and the same-index leading-k
                          span), its best-matching run-B PC and that |r| -
                          rotation shows as containment ~ 1 with a poor
                          matched-index |r|.
  pc_crosscorr.png        |correlation| heatmap, A x B - rotations appear as
                          small blocks around the diagonal, swaps as
                          off-diagonal dots.
  pc_set_concordance.md   the summary with interpretation.

Usage:
  python3 pc_set_concordance.py --normal svd.pcs.normal.txt --fast svd.pcs.fast.txt
"""

import argparse
import os

import numpy as np


def read_pcs(path):
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        samples, rows = [], []
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            samples.append(fields[0])
            rows.append([float(v) for v in fields[1:]])
    return samples, np.asarray(rows), len(header) - 1


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--normal", required=True, help="svd.pcs.txt from run A (the reference)")
    parser.add_argument("--fast", required=True, help="svd.pcs.txt from run B (the comparison)")
    parser.add_argument("--label-a", default="normal", help="display name for the --normal run")
    parser.add_argument("--label-b", default="fast", help="display name for the --fast run")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    label_a, label_b = args.label_a, args.label_b

    n_samples, n_mat, n_k = read_pcs(args.normal)
    f_samples, f_mat, f_k = read_pcs(args.fast)
    f_index = {s: i for i, s in enumerate(f_samples)}
    shared = [s for s in n_samples if s in f_index]
    A = n_mat[[n_samples.index(s) for s in shared], :] if shared != n_samples else n_mat
    B = f_mat[[f_index[s] for s in shared], :]
    k = min(n_k, f_k)
    A, B = A[:, :k], B[:, :k]
    print(f"{len(shared)} shared samples, {k} PCs per run")

    # Columns are right singular vectors and should be orthonormal already;
    # orthonormalize anyway so nothing below depends on that assumption.
    Qa, _ = np.linalg.qr(A)
    Qb, _ = np.linalg.qr(B)

    # ── Principal angles between nested leading-k subspaces ─────────────────
    ks = sorted({2, 5, 10, 20, 50, 100, 150, k} - {1})
    ks = [x for x in ks if x <= k]
    lines = ["k\tmin_cos\tlargest_angle_deg\tmean_cos\tn_cos_below_0.99\tn_cos_below_0.90"]
    angle_rows = []
    for kk in ks:
        cosines = np.linalg.svd(Qa[:, :kk].T @ Qb[:, :kk], compute_uv=False)
        cosines = np.clip(cosines, -1, 1)
        worst = float(cosines.min())
        angle = float(np.degrees(np.arccos(worst)))
        n99, n90 = int((cosines < 0.99).sum()), int((cosines < 0.90).sum())
        angle_rows.append((kk, worst, angle, float(cosines.mean()), n99, n90))
        lines.append(f"{kk}\t{worst:.6f}\t{angle:.3f}\t{cosines.mean():.6f}\t{n99}\t{n90}")
    with open(os.path.join(args.out_dir, "subspace_angles.tsv"), "w") as out:
        out.write("\n".join(lines) + "\n")

    # ── Cross-correlation and per-PC containment ─────────────────────────────
    Az = (A - A.mean(0)) / A.std(0)
    Bz = (B - B.mean(0)) / B.std(0)
    R = (Az.T @ Bz) / len(shared)                       # k x k correlations
    absR = np.abs(R)
    best_match = absR.argmax(axis=1)
    best_r = absR.max(axis=1)
    swaps = int((best_match != np.arange(k)).sum())

    # Containment: squared projection of each normal PC onto fast spans
    proj_full = Qb.T @ Qa                                # coordinates in fast basis
    contain_full = (proj_full ** 2).sum(axis=0)          # onto all k fast PCs
    contain_lead = np.array([
        (proj_full[: i + 1, i] ** 2).sum() for i in range(k)   # onto fast PC1..i+1
    ])
    with open(os.path.join(args.out_dir, "pc_containment.tsv"), "w") as out:
        out.write("PC\tmatched_abs_r\tbest_match_other_PC\tbest_abs_r\t"
                  "containment_in_all_other\tcontainment_in_leading_same_k\n")
        for i in range(k):
            out.write(f"{i+1}\t{absR[i, i]:.6f}\t{best_match[i]+1}\t{best_r[i]:.6f}\t"
                      f"{contain_full[i]:.6f}\t{contain_lead[i]:.6f}\n")

    # ── Heatmap ──────────────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8.5, 7.5), constrained_layout=True)
        im = ax.imshow(absR, cmap="viridis", vmin=0, vmax=1, origin="upper",
                       interpolation="nearest")
        ax.set_xlabel(f"{label_b} PC")
        ax.set_ylabel(f"{label_a} PC")
        ax.set_title(f"|correlation|, {label_a} x {label_b} - rotation blocks hug the diagonal")
        fig.colorbar(im, ax=ax, shrink=0.85)
        fig.savefig(os.path.join(args.out_dir, "pc_crosscorr.png"), dpi=150)
    except ImportError:
        print("matplotlib not available - skipping pc_crosscorr.png")

    # ── Summary ──────────────────────────────────────────────────────────────
    diag = np.diag(absR)
    md = [f"# PC set concordance, {label_a} vs {label_b}", ""]
    md.append(f"{len(shared)} shared samples; {k} PCs per run.")
    md.append("")
    md.append("| leading k | largest principal angle | min canonical correlation | cosines < 0.99 |")
    md.append("|---|---|---|---|")
    for kk, worst, angle, _, n99, _ in angle_rows:
        md.append(f"| {kk} | {angle:.2f}° | {worst:.6f} | {n99} of {kk} |")
    md.append("")
    md.append(f"- Matched-index |r|: median {np.median(diag):.4f}, "
              f"{int((diag < 0.9).sum())} of {k} below 0.9.")
    md.append(f"- Best-match |r|: median {np.median(best_r):.4f}, "
              f"{int((best_r < 0.9).sum())} below 0.9; {swaps} PCs match a different index.")
    md.append(f"- Containment of each {label_a} PC in the span of all {k} {label_b} PCs: "
              f"median {np.median(contain_full):.4f}, minimum {contain_full.min():.4f}.")
    md.append("")
    md.append("Reading it: matched-index |r| asks 'is PC i the same vector', containment asks "
              "'does PC i live inside the other run's space', and the principal angles ask "
              "'could any downstream analysis that uses the leading k as a set - regression "
              "covariates included - tell the runs apart'. Rotations within near-degenerate "
              "clusters destroy the first, and leave the last two untouched.")
    with open(os.path.join(args.out_dir, "pc_set_concordance.md"), "w") as out:
        out.write("\n".join(md) + "\n")

    print("\n".join(md[4:4 + len(ks) + 1]))
    print(f"\nwrote subspace_angles.tsv, pc_containment.tsv, pc_set_concordance.md"
          f"{', pc_crosscorr.png' if os.path.exists(os.path.join(args.out_dir, 'pc_crosscorr.png')) else ''}")


if __name__ == "__main__":
    main()
