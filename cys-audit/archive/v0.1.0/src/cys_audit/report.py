"""Writers: audit.json, summary.csv, audit.html, figures/*.svg, manifest.json.

Outputs contain no timestamps and no absolute paths, so two runs on the same input with the same
seed produce byte-identical files (the manifest records sha256 for each).
"""
from __future__ import annotations

import csv
import html
import json
import os

from .audit import sha256_file
from .checks import ORDER

COLOURS = {"PASS": "#2e7d32", "WARNING": "#b26a00", "FAIL": "#c62828", "UNDECIDABLE": "#616161",
           "ATTENUATED": "#b26a00", "NOT_EVALUATED": "#616161"}


def write_json(rec, path):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(rec, fh, indent=1, ensure_ascii=False, allow_nan=False, sort_keys=False)
        fh.write("\n")


def summary_rows(rec):
    rows = []
    for name in ORDER:
        t = rec["tests"].get(name)
        if t is None:
            continue
        ci = t.get("ci") or [None, None]
        rows.append({"dataset": rec["dataset"], "test": name, "status": t["status"], "estimate": t.get("estimate"),
                     "ci_low": ci[0], "ci_high": ci[1], "null": t.get("null"),
                     "margin": t.get("margin") if not isinstance(t.get("margin"), dict) else "see audit.json",
                     "p_value": t.get("p_value"), "q_value_bh": t.get("q_value_bh"),
                     "n_positive": t.get("n_positive"), "n_background": t.get("n_background"),
                     "reason": t.get("reason")})
    return rows


def write_csv(rec, path):
    rows = summary_rows(rec)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def interval_svg(t, width=420, height=70):
    """Point + interval against the null and the equivalence margin, for interval-based tests."""
    if t.get("ci") is None or t.get("estimate") is None or t.get("null") is None or not isinstance(t.get("margin"), (int, float)):
        return None
    lo, hi = t["ci"]
    if lo is None or hi is None:
        return None
    null, m, est = t["null"], t["margin"], t["estimate"]
    vmin = min(lo, null - 1.5 * m, est)
    vmax = max(hi, null + 1.5 * m, est)
    pad = 0.08 * (vmax - vmin or 1)
    vmin, vmax = vmin - pad, vmax + pad
    X = lambda v: 20 + (v - vmin) / (vmax - vmin) * (width - 40)
    col = COLOURS.get(t["status"], "#000")
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
             f'font-family="Arial, Helvetica, sans-serif" font-size="11">',
             f'<rect x="{X(null - m):.1f}" y="10" width="{X(null + m) - X(null - m):.1f}" height="30" fill="#e8f5e9"/>',
             f'<line x1="{X(null):.1f}" y1="6" x2="{X(null):.1f}" y2="44" stroke="#555" stroke-dasharray="3,2"/>',
             f'<line x1="{X(lo):.1f}" y1="25" x2="{X(hi):.1f}" y2="25" stroke="{col}" stroke-width="2.5"/>',
             f'<circle cx="{X(est):.1f}" cy="25" r="4.5" fill="{col}"/>',
             f'<text x="{X(null):.1f}" y="60" text-anchor="middle" fill="#555">null {null:g}</text>',
             f'<text x="20" y="60" fill="#555">shaded: +/- margin {m:g}</text>',
             "</svg>"]
    return "\n".join(parts)


def write_figures(rec, fig_dir):
    os.makedirs(fig_dir, exist_ok=True)
    written = []
    for name in ORDER:
        t = rec["tests"].get(name)
        if t is None:
            continue
        svg = interval_svg(t)
        if svg:
            p = os.path.join(fig_dir, f"{name}.svg")
            with open(p, "w", encoding="utf-8") as fh:
                fh.write(svg + "\n")
            written.append(p)
    return written


def write_html(rec, path, fig_dir):
    e = html.escape
    rows = []
    for name in ORDER:
        t = rec["tests"].get(name)
        if t is None:
            continue
        ci = t.get("ci")
        ci_s = "" if not ci or ci[0] is None else f"[{ci[0]:.4g}, {ci[1]:.4g}]"
        est = "" if t.get("estimate") is None else f"{t['estimate']:.4g}"
        svg_p = os.path.join(fig_dir, f"{name}.svg")
        svg = open(svg_p, encoding="utf-8").read() if os.path.exists(svg_p) else ""
        rows.append(f"<tr><td><b>{e(name)}</b></td><td style='color:{COLOURS[t['status']]};font-weight:bold'>"
                    f"{e(t['status'])}</td><td>{e(est)}</td><td>{e(ci_s)}</td><td>{e(str(t.get('n_positive')))}"
                    f" / {e(str(t.get('n_background')))}</td><td>{e(t.get('reason') or '')}<br>{svg}</td></tr>")
    claim = rec["claim"]
    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>Cys-Audit report - {e(rec['dataset'])}</title>
<style>body{{font-family:Arial,Helvetica,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px;color:#222}}
table{{border-collapse:collapse;width:100%}}td,th{{border-bottom:1px solid #ddd;padding:8px;vertical-align:top;text-align:left}}
code{{background:#f4f4f4;padding:1px 4px}}</style></head><body>
<h1>Cys-Audit report</h1>
<p>Dataset <code>{e(rec['dataset'])}</code> &middot; modification <code>{e(str(rec['modification']))}</code> &middot;
protease <code>{e(rec['protease'])}</code> &middot; background <code>{e(rec['background'])}</code> &middot;
tool {e(rec['version'])} &middot; seed {rec['seed']}</p>
<p>Rows {rec['input']['n_rows']}, positives {rec['input']['n_positive']}, detected {rec['input']['n_detected']},
proteins {rec['input']['n_proteins']}, clusters {rec['input']['n_clusters']}.</p>
<h2>Overall claim status: <span style="color:{COLOURS.get(rec['overall_claim_status'], '#000')}">{e(rec['overall_claim_status'])}</span></h2>
<p>{e(claim.get('reason', ''))}</p>
<h2>Artefact tests</h2>
<table><tr><th>Test</th><th>Status</th><th>Estimate</th><th>Interval</th><th>n pos / bg</th><th>Reason</th></tr>
{''.join(rows)}</table>
<h2>How to read a status</h2>
<ul><li><b>PASS</b>: the interval lies inside the equivalence margin - an artefact of meaningful size is excluded.</li>
<li><b>WARNING</b>: the interval excludes the null but does not establish a meaningful size.</li>
<li><b>FAIL</b>: the interval lies beyond the margin - an artefact of meaningful size is established.</li>
<li><b>UNDECIDABLE</b>: inputs missing, classes too small, or the interval is too wide to decide.</li></ul>
<p>A status describes the dataset on one axis. It is not a verdict on any biological claim unless a claim
feature was supplied (overall claim status above). Notes: {e('; '.join(rec['input']['notes']) or 'none')}</p>
</body></html>
"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(doc)


def write_all(rec, out_dir, input_paths, command):
    os.makedirs(out_dir, exist_ok=True)
    fig_dir = os.path.join(out_dir, "figures")
    write_json(rec, os.path.join(out_dir, "audit.json"))
    write_csv(rec, os.path.join(out_dir, "summary.csv"))
    figs = write_figures(rec, fig_dir)
    write_html(rec, os.path.join(out_dir, "audit.html"), fig_dir)
    outputs = ["audit.json", "summary.csv", "audit.html"] + [os.path.relpath(f, out_dir) for f in figs]
    manifest = {"tool": "cys-audit", "version": rec["version"], "command": command,
                "inputs": {os.path.basename(p): sha256_file(p) for p in input_paths if p},
                "outputs": {o: sha256_file(os.path.join(out_dir, o)) for o in sorted(outputs)},
                "environment": rec["environment"]}
    write_json(manifest, os.path.join(out_dir, "manifest.json"))
    return manifest
