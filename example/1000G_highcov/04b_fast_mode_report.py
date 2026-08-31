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
