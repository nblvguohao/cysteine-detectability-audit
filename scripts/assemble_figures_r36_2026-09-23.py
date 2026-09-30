"""Assemble figures_r36/ and source_data_r36/ for R36 (standard library only).

R36 = R35 with Figure 1's panel a redrawn to drop the per-step artefact number (scripts/plot_fig1_biorender_r36_2026-09-23.py,
already written into figures_r36/). Figs 2-7 and their Source Data are byte copies of figures_r35/ and source_data_r35/
(sha256 asserted) -- nothing about them changed this round. Figure 1's Source Data table is regenerated from this
round's own audit (results/fig1_biorender_r36_2026-09-23_audit.json printed_values); the printed counts are identical
to R35's (same source tables, same assertions), so the table's values are unchanged, only its provenance pointer is.
Refuses if source_data_r36/ already exists or if any figures_r36 Fig 2-7 file already exists.
"""
import csv, hashlib, json, pathlib, shutil, sys
ROOT = pathlib.Path(__file__).resolve().parent.parent
sha = lambda p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
fig_src, fig_dst = ROOT / "figures_r35", ROOT / "figures_r36"
sd_src, sd_dst = ROOT / "source_data_r35", ROOT / "source_data_r36"
if sd_dst.exists() or any(fig_dst.glob("Fig[2-7]_*")):
    sys.exit("REFUSE: outputs exist")
if not (fig_dst / "Fig1_overview.pdf").exists():
    sys.exit("REFUSE: run plot_fig1_biorender_r36 first")
copied = {}
for f in sorted(fig_src.glob("Fig[2-7]_*")):
    shutil.copy2(f, fig_dst / f.name); assert sha(f) == sha(fig_dst / f.name); copied[str((fig_dst / f.name).relative_to(ROOT))] = sha(f)
sd_dst.mkdir()
for f in sorted(sd_src.iterdir()):
    if f.is_file() and f.name != "Source_Data_Fig1_overview.csv":
        shutil.copy2(f, sd_dst / f.name); assert sha(f) == sha(sd_dst / f.name); copied[str((sd_dst / f.name).relative_to(ROOT))] = sha(f)
aud = json.loads((ROOT / "results/fig1_biorender_r36_2026-09-23_audit.json").read_text())
assert aud["all_pass"]
pv = aud["printed_values"]
rows = [dict(figure="Fig1", panel="b", row=l, field="records", value=n, source_table=aud["sources"]["survey"]) for l, n in pv["funnel"]]
rows += [dict(figure="Fig1", panel="c", row=l.replace("\n", " "), field="classified_analyses", value=n, source_table=aud["sources"]["survey"]) for l, n in pv["composition"]]
rows += [dict(figure="Fig1", panel="d", row=l, field="claims", value=n, source_table=aud["sources"]["tally"]) for l, n in pv["verdicts"]]
out = sd_dst / "Source_Data_Fig1_overview.csv"
with out.open("w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
# cross-check: R36's Fig1 values equal R35's Fig1 values (only panel a's drawing changed, not the data)
prior = list(csv.DictReader((ROOT / "source_data_r35/Source_Data_Fig1_overview.csv").open(encoding="utf-8", newline="")))
current = list(csv.DictReader(out.open(encoding="utf-8", newline="")))
prior_vals = [(r["panel"], r["row"], r["field"], r["value"]) for r in prior]
current_vals = [(r["panel"], r["row"], r["field"], r["value"]) for r in current]
assert prior_vals == current_vals, "Fig1 Source Data values changed between R35 and R36; only the drawing should differ"
(ROOT / "results/figures_r36_assembly_2026-09-23_audit.json").write_text(json.dumps(dict(
    script=pathlib.Path(__file__).name, script_sha256=sha(__file__), byte_copies=copied,
    fig1_source_data={str(out.relative_to(ROOT)): sha(out)}, fig1_rows=len(rows),
    fig1_values_unchanged_from_r35=True,
    fig1_outputs={k: v for k, v in aud["outputs"].items()}), indent=1))
print("copied", len(copied), "fig1 source rows", len(rows), "values_unchanged_from_r35", True)
