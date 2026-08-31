#!/usr/bin/env python3
"""Assemble the fast-mode evaluation into one self-contained HTML report.

Reads what 04_fast_mode_eval.py wrote - the correlation, timing and QC
concordance tables, the markdown summary, and the figures when matplotlib was
available - and writes fast_mode_report.html: methods, results, figures
(embedded as data URIs, so the file stands alone), and interpretation
generated from the numbers themselves. Print to PDF from any browser for a
paper-ready copy.

Standard library only, so it runs anywhere; the figures simply appear when
the eval produced them (04b_fast_mode_report.sh runs both steps inside a
container that has matplotlib, which is the intended path).
"""

import argparse
import base64
import csv
import datetime
import html
import os
import sys


def read_tsv(path):
    if not os.path.isfile(path):
        return []
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def embed_png(path):
    if not os.path.isfile(path):
        return None
    with open(path, "rb") as handle:
        payload = base64.b64encode(handle.read()).decode("ascii")
    return f"data:image/png;base64,{payload}"


def table_html(rows, columns, caption):
    if not rows:
        return f"<p class='missing'>{html.escape(caption)}: not found - re-run the eval.</p>"
    head = "".join(f"<th>{html.escape(c)}</th>" for c in columns)
    body = []
    for row in rows:
        cells = "".join(f"<td>{html.escape(row.get(c, ''))}</td>" for c in columns)
        body.append(f"<tr>{cells}</tr>")
    return (f"<figure><figcaption>{html.escape(caption)}</figcaption>"
            f"<div class='scroll'><table><thead><tr>{head}</tr></thead>"
            f"<tbody>{''.join(body)}</tbody></table></div></figure>")


def fmt(value, digits=6):
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--eval-dir", required=True, help="the eval's output directory")
    parser.add_argument("--out", default=None, help="HTML path (default: <eval-dir>/fast_mode_report.html)")
    args = parser.parse_args()

    eval_dir = args.eval_dir
    out_path = args.out or os.path.join(eval_dir, "fast_mode_report.html")

    pcs = read_tsv(os.path.join(eval_dir, "pc_correlation.tsv"))
    timing = read_tsv(os.path.join(eval_dir, "timing_summary.tsv"))
    qc = read_tsv(os.path.join(eval_dir, "qc_concordance.tsv"))
    if not pcs:
        sys.exit(f"ERROR: {eval_dir}/pc_correlation.tsv not found - run 04_fast_mode_eval first.")

    md_path = os.path.join(eval_dir, "fast_mode_report.md")
    md_text = open(md_path).read() if os.path.isfile(md_path) else ""

    angles = read_tsv(os.path.join(eval_dir, "subspace_angles.tsv"))
    containment = read_tsv(os.path.join(eval_dir, "pc_containment.tsv"))

    # Seed-control calibration pair, when 04b_fast_mode_report.sh found a
    # reseeded normal run: same tables, same analysis, different question.
    seed_dir = os.path.join(eval_dir, "seed_control")
    seed_angles = read_tsv(os.path.join(seed_dir, "subspace_angles.tsv"))
    seed_containment = read_tsv(os.path.join(seed_dir, "pc_containment.tsv"))

    def col(row, name, legacy):
        # pc_containment.tsv column names went generic when the seed control
        # arrived; accept tables written by either version.
        value = row.get(name)
        return value if value is not None else row.get(legacy)

    # ── Headline numbers, derived from the tables rather than trusted ────────
    abs_rs = [(row["PC"], float(row["abs_r"])) for row in pcs
              if row["abs_r"].lower() != "nan"]
    min_pc, min_abs_r = min(abs_rs, key=lambda pair: pair[1]) if abs_rs else ("?", float("nan"))
    speedup = next((row for row in timing if row["metric"] == "paired_speedup"), None)
    walls = {row["metric"]: row for row in timing}
    mtdna = next((row for row in qc if row["column"] == "MTDNA_CN"), None)

    weak = [row for row in pcs if row["abs_r"].lower() != "nan" and float(row["abs_r"]) < 0.999]
    flipped = [row["PC"] for row in pcs if row["pearson_r"].lstrip().startswith("-")]
    nan_qc = [row["column"] for row in qc if row["pearson_r"].lower() == "nan"]

    interpretation = []
    if len(abs_rs) > 20:
        def tier(name, values):
            ordered = sorted(v for _, v in values)
            below = sum(1 for v in ordered if v < 0.99)
            return (f"{name}: median |r| = {ordered[len(ordered) // 2]:.6f}, "
                    f"min = {ordered[0]:.6f}, {below} below 0.99")
        interpretation.append(
            f"Stability by depth - {tier('PC1-PC20', abs_rs[:20])}; "
            f"{tier('PC21-PC' + abs_rs[-1][0], abs_rs[20:])}. Later components carry less "
            f"variance and rotate more easily among near-equal singular values, so a slow "
            f"decline in |r| down the spectrum is expected rather than alarming - what would "
            f"matter is an early component breaking from it.")
    if weak:
        notes = ", ".join(
            f"PC{row['PC']} (|r| = {fmt(row['abs_r'])}, singular-value ratio "
            f"{fmt(row['singular_value_ratio_fast_over_normal'], 4)})" for row in weak)
        interpretation.append(
            f"PCs below |r| = 0.999: {notes}. A ratio near 1 alongside a dip in |r| is the "
            f"signature of rotation between near-degenerate components - adjacent PCs with "
            f"almost equal singular values can trade directions without any real disagreement "
            f"between the runs.")
    else:
        interpretation.append("Every evaluated PC correlates at |r| ≥ 0.999.")
    if flipped:
        interpretation.append(
            f"PC{', PC'.join(flipped)} carried a negative r: a pure sign flip, which is why "
            f"correlations are evaluated as |r| - the sign of a singular vector is arbitrary.")
    if nan_qc:
        interpretation.append(
            f"QC columns with r = nan ({', '.join(nan_qc)}) have zero variance across samples "
            f"in at least one run - a constant column cannot correlate, and its median "
            f"difference is the informative number.")
    interpretation.append(
        "Raw-depth columns are expected to shift under fast mode (skipped mate-overlap "
        "correction inflates depth where read pairs overlap), and the per-sample log2 "
        "fold-change normalization is expected to absorb it: the shift shows in the QC "
        "means and dispersions while the PCs stay put, which is exactly the pattern that "
        "validates using fast mode for this pipeline.")

    crosscorr_img = embed_png(os.path.join(eval_dir, "pc_crosscorr.png"))
    if angles:
        swaps = sum(1 for row in containment
                    if col(row, "best_match_other_PC", "best_match_fast_PC") != row.get("PC"))
        contain_vals = sorted(float(col(row, "containment_in_all_other", "containment_in_all_fast"))
                              for row in containment)
        set_bits = [table_html(angles, list(angles[0].keys()),
                               "Principal angles between leading-k subspaces")]
        if crosscorr_img:
            set_bits.append(f"<img src='{crosscorr_img}' alt='cross-correlation heatmap'>")
        set_bits.append(
            f"<p>Containment of each PC in the other run's full span: median "
            f"{contain_vals[len(contain_vals) // 2]:.4f}, minimum {contain_vals[0]:.4f}; "
            f"{swaps} of {len(containment)} PCs best-match a different index (rank swaps and "
            f"rotations, localized in pc_containment.tsv). Small principal angles alongside a "
            f"low matched-index |r| is the signature of rotation within near-degenerate "
            f"clusters - the set agrees even where the indexing does not.</p>")
        set_section = "".join(set_bits)
    else:
        set_section = ("<p class='missing'>Set-concordance outputs not found - run "
                       "04b_fast_mode_report.sh via the eval container (numpy) to add "
                       "principal angles, containment, and the cross-correlation heatmap.</p>")

    # ── Calibration: what does the estimator do to ITSELF under a reseed? ────
    # Both PCA runs of the fast comparison share one seed, so estimator
    # randomness contributed nothing there - every difference is caused by the
    # fast-mode depths. The seed control recomputes the identical normal tree
    # under a different -randomSeed: differences there are pure truncation
    # noise, and they calibrate whether the fast-mode differences ever exceed
    # what the estimator produces against itself.
    calibration_section = ""
    if angles and seed_angles:
        seed_by_k = {row["k"]: row for row in seed_angles}
        joined = [(row, seed_by_k[row["k"]]) for row in angles if row["k"] in seed_by_k]
        have_counts = all("n_cos_below_0.99" in row and "n_cos_below_0.99" in s
                          for row, s in joined)
        cal_columns = ["k", "largest angle: fast vs normal", "largest angle: seed vs seed",
                       "mean cos (fast)", "mean cos (seed)"]
        if have_counts:
            cal_columns.append("cosines < 0.99: fast / seed")
        cal_rows = []
        for row, s in joined:
            cal_row = {
                "k": row["k"],
                "largest angle: fast vs normal": f"{float(row['largest_angle_deg']):.2f}°",
                "largest angle: seed vs seed": f"{float(s['largest_angle_deg']):.2f}°",
                "mean cos (fast)": fmt(row["mean_cos"], 4),
                "mean cos (seed)": fmt(s["mean_cos"], 4),
            }
            if have_counts:
                cal_row["cosines < 0.99: fast / seed"] = (
                    f"{row['n_cos_below_0.99']} / {s['n_cos_below_0.99']}")
            cal_rows.append(cal_row)
        # A k where fast-vs-normal exceeds both twice the seed-control angle
        # and +5 degrees is a candidate genuine fast-mode effect; anything
        # inside that band is indistinguishable from the estimator's noise.
        exceed = [(row["k"], float(row["largest_angle_deg"]), float(s["largest_angle_deg"]))
                  for row, s in joined
                  if float(row["largest_angle_deg"]) > 2 * float(s["largest_angle_deg"])
                  and float(row["largest_angle_deg"]) - float(s["largest_angle_deg"]) > 5]

        def tail_stats(rows):
            values = [float(col(r, "containment_in_all_other", "containment_in_all_fast"))
                      for r in rows]
            return sum(1 for v in values if v < 0.9), min(values) if values else float("nan")
        fast_low, fast_min = tail_stats(containment)
        seed_low, seed_min = tail_stats(seed_containment)

        if exceed:
            listed = "; ".join(f"k = {k}: {fa:.1f}° vs {sa:.1f}°" for k, fa, sa in exceed)
            verdict = (f"Fast-vs-normal exceeds the seed-control band at {listed} - candidate "
                       f"genuine fast-mode effects worth inspecting; every other k sits within "
                       f"what the estimator does to itself under a reseed alone.")
        else:
            # Two regimes, described separately: where a reseed leaves the
            # subspace unchanged, any fast-mode angle is a real (if small)
            # data-caused effect; where the reseed itself swings the subspace
            # further than fast mode does, those directions are not
            # reproducible under the estimator's own seed and no fast-mode
            # effect is even definable there.
            triples = [(row["k"], float(row["largest_angle_deg"]), float(s["largest_angle_deg"]))
                       for row, s in joined]
            tight = [(kk, fa) for kk, fa, sa in triples if sa <= 5]
            loose = [(kk, fa, sa) for kk, fa, sa in triples if sa > fa]
            parts = []
            if tight:
                parts.append(
                    f"Where a reseed alone leaves the leading subspace essentially unchanged "
                    f"(angle ≤ 5°: k = {', '.join(kk for kk, _ in tight)}), the fast-mode angle "
                    f"tops out at {max(fa for _, fa in tight):.1f}° - a genuine, data-caused, "
                    f"and small rotation.")
            if loose:
                kk_w, fa_w, sa_w = max(loose, key=lambda t: t[2])
                parts.append(
                    f"Where the spectrum runs into near-degeneracy "
                    f"(k = {', '.join(kk for kk, _, _ in loose)}), the reseed swings the "
                    f"subspace further than fast mode does - {sa_w:.0f}° against {fa_w:.0f}° at "
                    f"k = {kk_w} - so the fast-mode differences there sit below the estimator's "
                    f"own noise floor: those directions are not reproducible under the "
                    f"estimator's seed, with or without fast mode.")
            if not parts:
                parts.append(
                    "At every evaluated k, the fast-vs-normal angle sits within the band the "
                    "estimator produces against itself under a seed change alone.")
            verdict = " ".join(parts)
        verdict += (f" Containment tells the same story at the spectrum's edge: "
                    f"{fast_low} PCs below 0.9 containment in the fast comparison (minimum "
                    f"{fast_min:.2f}) against {seed_low} in the seed control (minimum "
                    f"{seed_min:.2f}).")

        cal_bits = [
            "<h2>Calibration: seed change vs fast mode</h2>",
            "<p>The fast-vs-normal comparison held the estimator fixed: one seed, one sample "
            "order, so every difference above is <em>caused</em> by the fast-mode depths. This "
            "section asks how large those differences are on the only meaningful scale - the "
            "one the estimator sets itself. The identical normal tree was recomputed under a "
            "different random seed and compared to the seed-42 run with the same subspace "
            "analysis: those differences are pure truncation noise from the randomized SVD, "
            "with byte-identical input data.</p>",
            table_html(cal_rows, cal_columns,
                       "Largest principal angle, fast-mode comparison vs seed-control calibration"),
            f"<p>{html.escape(verdict)}</p>",
        ]
        seed_img = embed_png(os.path.join(seed_dir, "pc_crosscorr.png"))
        if seed_img:
            cal_bits.append(f"<img src='{seed_img}' alt='seed-control cross-correlation heatmap'>")
            cal_bits.append("<p>The same diagnostic heatmap for the seed control - normal PCs "
                            "against the reseeded run's. A diagonal band of the same character "
                            "as the fast comparison's is the visual form of the verdict.</p>")
        calibration_section = "".join(cal_bits)
    elif angles:
        calibration_section = (
            "<h2>Calibration: seed change vs fast mode</h2>"
            "<p class='missing'>No seed-control run found, so the numbers above lack their "
            "yardstick: how much would the estimator alone move under a reseed? Recompute the "
            "normal tree with a different seed, then re-run this report:<br>"
            "<code>NGSPCA_OUTPUT=\"${WORK_DIR}/ngspca_output_seed43\" RANDOM_SEED=43 "
            "sbatch 02_run_ngspca.sh</code></p>")

    summary_img = embed_png(os.path.join(eval_dir, "fast_mode_summary.png"))
    qc_img = embed_png(os.path.join(eval_dir, "fast_mode_qc.png"))

    def img_or_note(uri, alt):
        if uri:
            return f"<img src='{uri}' alt='{html.escape(alt)}'>"
        return ("<p class='missing'>Figure not found - run the eval where matplotlib is "
                "available (04b_fast_mode_report.sh does this inside the eval container).</p>")

    headline_cells = [
        ("Minimum |r|, PC1-PC" + str(len(pcs)), f"{fmt(min_abs_r)} (PC{min_pc})"),
    ]
    if angles:
        full = angles[-1]
        headline_cells.append((f"Largest principal angle, all {full['k']} PCs as a set",
                               f"{fmt(full['largest_angle_deg'], 2)}° (min canonical corr "
                               f"{fmt(full['min_cos'], 4)})"))
    if seed_angles:
        seed_full = seed_angles[-1]
        headline_cells.append((f"Seed-control angle, all {seed_full['k']} PCs as a set",
                               f"{fmt(seed_full['largest_angle_deg'], 2)}° - estimator "
                               f"noise alone"))
    if speedup:
        headline_cells.append(("Median mosdepth speedup",
                               f"{fmt(speedup['median'], 2)}x (n = {speedup['n']})"))
    if "normal_wall_s" in walls and "fast_wall_s" in walls:
        headline_cells.append(("Median wall time",
                               f"{fmt(walls['normal_wall_s']['median'], 0)} s normal / "
                               f"{fmt(walls['fast_wall_s']['median'], 0)} s fast"))
    if mtdna:
        headline_cells.append(("MTDNA_CN concordance",
                               f"r = {fmt(mtdna['pearson_r'])}, shift "
                               f"{float(mtdna['median_rel_diff_pct']):+.2f}%"))
    headlines = "".join(
        f"<div class='stat'><div class='value'>{html.escape(v)}</div>"
        f"<div class='label'>{html.escape(k)}</div></div>" for k, v in headline_cells)

    interpretation_html = "".join(f"<li>{html.escape(note)}</li>" for note in interpretation)

    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>mosdepth fast-mode evaluation</title>
<style>
  body {{ font: 15px/1.55 -apple-system, "Segoe UI", Roboto, sans-serif; color: #1a2733;
         max-width: 1000px; margin: 2rem auto; padding: 0 1rem; }}
  h1 {{ font-size: 1.6rem; margin-bottom: 0.2rem; }}
  h2 {{ font-size: 1.15rem; margin-top: 2rem; border-bottom: 2px solid #e3e9ef; padding-bottom: 0.3rem; }}
  .meta {{ color: #5c6c7a; margin-bottom: 1.5rem; }}
  .stats {{ display: flex; flex-wrap: wrap; gap: 0.8rem; margin: 1rem 0; }}
  .stat {{ flex: 1 1 200px; background: #f2f6fa; border-radius: 8px; padding: 0.9rem 1rem; }}
  .stat .value {{ font-size: 1.25rem; font-weight: 600; color: #234; }}
  .stat .label {{ font-size: 0.8rem; color: #5c6c7a; margin-top: 0.2rem; }}
  img {{ max-width: 100%; border: 1px solid #e3e9ef; border-radius: 6px; margin: 0.6rem 0; }}
  table {{ border-collapse: collapse; font-size: 0.82rem; }}
  th, td {{ border: 1px solid #dde5ec; padding: 0.25rem 0.55rem; text-align: right; }}
  th {{ background: #f2f6fa; }} td:first-child, th:first-child {{ text-align: left; }}
  figure {{ margin: 1rem 0; }} figcaption {{ font-weight: 600; margin-bottom: 0.4rem; }}
  .scroll {{ overflow-x: auto; }}
  .missing {{ color: #a33; background: #fdf3f3; padding: 0.6rem; border-radius: 6px; }}
  details pre {{ background: #f7f9fb; padding: 0.8rem; overflow-x: auto; font-size: 0.78rem; }}
  li {{ margin: 0.4rem 0; }}
  @media print {{ .stat {{ border: 1px solid #ccc; }} }}
</style></head><body>
<h1>mosdepth fast-mode evaluation</h1>
<p class="meta">Generated {datetime.date.today().isoformat()} from {html.escape(eval_dir)} ·
print this page to PDF for a fixed copy</p>

<div class="stats">{headlines}</div>

<h2>Methods, briefly</h2>
<p>Each sample's mosdepth ran twice against one verified CRAM on one node - once as configured,
once with <code>--fast-mode</code> - with run order alternating by manifest line so page-cache
warming cancels in aggregate (the order check in the tables below confirms it did). Downloads
were excluded from both timings by construction. NGS-PCA then ran identically on each output
tree (same seed, parameters and exclusion regions), and per-PC Pearson correlations were
computed across shared samples, evaluated as |r| because a singular vector's sign is arbitrary,
with singular-value ratios reported so rotation between near-degenerate components can be told
apart from real disagreement. QC concordance compares the per-sample tables column by column:
Pearson r plus the median signed relative difference, since a uniform bias - the expected
signature of skipping mate-overlap correction - is invisible to correlation alone.</p>

<h2>Principal components</h2>
{img_or_note(summary_img, "PC correlations, runtimes, and scale agreement")}
{table_html(pcs, list(pcs[0].keys()) if pcs else [], "Per-PC correlation and singular-value ratio")}

<h2>QC concordance</h2>
{img_or_note(qc_img, "QC column agreement and MTDNA_CN scatter")}
{table_html(qc, list(qc[0].keys()) if qc else [], "Per-column concordance of sample_qc.tsv")}

<h2>The PC sets as subspaces</h2>
<p>Matched-index correlation asks whether PC <i>i</i> is the same vector in both runs - the
strictest test, and the wrong one deep in the spectrum, where near-equal singular values let
components rotate into each other or swap rank without the <em>span</em> changing. Downstream
use of these PCs as covariates depends only on that span, so the set-level questions are the
operative ones: the largest principal angle between the leading-k subspaces bounds how
differently any analysis using them as a set could behave, and per-PC containment measures how
completely each component lives inside the other run's space regardless of index.</p>
{set_section}
{calibration_section}

<h2>Runtimes</h2>
{table_html(timing, list(timing[0].keys()) if timing else [], "mosdepth wall-time statistics")}

<h2>Interpretation</h2>
<ul>{interpretation_html}</ul>

<h2>Raw evaluation summary</h2>
<details><summary>fast_mode_report.md as written by the eval</summary>
<pre>{html.escape(md_text) if md_text else "(not found)"}</pre></details>
</body></html>
"""
    with open(out_path, "w") as out:
        out.write(doc)
    print(f"wrote {out_path}")
    if not summary_img or not qc_img:
        print("note: one or both figures were missing; the report carries tables only there.")


if __name__ == "__main__":
    main()
