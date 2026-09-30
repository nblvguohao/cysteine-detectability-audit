"""Assemble figures_r35/ and source_data_r35/ for R35 (standard library only).

R35 = R34 with a new Figure 1 (scripts/plot_fig1_biorender_r35_2026-09-23.py, already written into figures_r35/).
Figs 2-7 and their Source Data are byte copies of figures_r34/ and source_data_r34/ (sha256 asserted). Figure 1 gets
a Source Data table for the first time: every count printed in panels b-d, taken from the figure script's own audit
(results/fig1_biorender_r35_2026-09-23_audit.json printed_values), which recomputed them from the stored tables.
Refuses if source_data_r35/ already exists or if any figures_r35 Fig 2-7 file already exists.
"""
import csv, hashlib, json, pathlib, shutil, sys
ROOT = pathlib.Path(__file__).resolve().parent.parent
sha = lambda p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
fig_src, fig_dst = ROOT / "figures_r34", ROOT / "figures_r35"
sd_src, sd_dst = ROOT / "source_data_r34", ROOT / "source_data_r35"
if sd_dst.exists() or any(fig_dst.glob("Fig[2-7]_*")):
    sys.exit("REFUSE: outputs exist")
if not (fig_dst / "Fig1_overview.pdf").exists():
    sys.exit("REFUSE: run plot_fig1_biorender_r35 first")
copied = {}
for f in sorted(fig_src.glob("Fig[2-7]_*")):
    shutil.copy2(f, fig_dst / f.name); assert sha(f) == sha(fig_dst / f.name); copied[str((fig_dst / f.name).relative_to(ROOT))] = sha(f)
sd_dst.mkdir()
for f in sorted(sd_src.iterdir()):
    if f.is_file():
        shutil.copy2(f, sd_dst / f.name); assert sha(f) == sha(sd_dst / f.name); copied[str((sd_dst / f.name).relative_to(ROOT))] = sha(f)
aud = json.loads((ROOT / "results/fig1_biorender_r35_2026-09-23_audit.json").read_text())
assert aud["all_pass"]
pv = aud["printed_values"]
rows = [dict(figure="Fig1", panel="b", row=l, field="records", value=n, source_table=aud["sources"]["survey"]) for l, n in pv["funnel"]]
rows += [dict(figure="Fig1", panel="c", row=l.replace("\n", " "), field="classified_analyses", value=n, source_table=aud["sources"]["survey"]) for l, n in pv["composition"]]
rows += [dict(figure="Fig1", panel="d", row=l, field="claims", value=n, source_table=aud["sources"]["tally"]) for l, n in pv["verdicts"]]
out = sd_dst / "Source_Data_Fig1_overview.csv"
with out.open("w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
(ROOT / "results/figures_r35_assembly_2026-09-23_audit.json").write_text(json.dumps(dict(
    script=pathlib.Path(__file__).name, script_sha256=sha(__file__), byte_copies=copied,
    fig1_source_data={str(out.relative_to(ROOT)): sha(out)}, fig1_rows=len(rows),
    fig1_outputs={k: v for k, v in aud["outputs"].items()}), indent=1))
print("copied", len(copied), "fig1 source rows", len(rows))
