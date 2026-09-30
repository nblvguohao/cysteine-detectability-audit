"""Figure 1 for R36: R35's figure with the per-step artefact number dropped from panel a.

User request (2026-09-23): "沿流程读作 4、1、3、5、2 ... 能不能改成 1、2、3、4、5，这样更好看一些". Renumbering the
whole manuscript to match workflow order was offered as the alternative and declined (it would touch 17+ "Artefact N"
references, five Results headings, the Fig 4 colour map, every figure's numbering, Source Data file names and the
cover letter -- a change on the scale of R31). The adopted fix is narrower: panel a no longer prints a number at all,
so there is nothing that can read out of Results order. Each artefact keeps the identity it has everywhere else in
the paper through its colour alone (S.ARTEFACT[num], the same map Fig 4 and Fig 6 use) -- the artefact *name* is now
set in that colour, bold, where the numbered badge used to sit. Nothing else about panel a changes: same illustration
crop, same step names, same control wording, same ink span positions. Panels b-d are untouched (same source tables,
same recomputed-and-asserted counts as R35).

  a  ILLUSTRATION + text, as in R35, minus the numbered badge. The FancyBboxPatch badge and its digit are removed;
     the artefact name is coloured instead of set in INK, and left-aligned with the step's own ink (no more offset
     for the badge width).
  b  survey funnel     \
  c  background used    > unchanged from R35: same source tables, same asserted values, same axes.
  d  claim outcomes      /

GATES (fixed before the first run)
  G1 printed counts recomputed and equal to the manuscript's      G2 minimum line width >= 1.0 pt
  G3 font sizes {7, 8}                                             G4 zero text collisions, planted positive control
  G5 no text off canvas                                            G6 the illustration file matches PROVENANCE sha256
  G7 scipilot-figure-skill visual_qa: no FAIL; WARNs recorded      G8 the saved PDF embeds only Arial (checked after)
  G9 no digit 1-5 is drawn inside the ax_note annotation band (the thing this round removes)
MODES: --draft DIR writes only to DIR. Final writes figures_r36/Fig1_overview.{pdf,png,svg} and refuses if they exist.
INTERPRETER: /Users/lyuguohao/Documents/lnrna复现/.venv/bin/python
"""
import csv
import hashlib
import importlib.util
import json
import pathlib
import re
import sys
from collections import Counter

ROOT = pathlib.Path(__file__).resolve().parent.parent
ILLU = ROOT / "external/biorender_fig1_20260923/fig1_workflow_biorender_v1.jpg"
PROV = ROOT / "external/biorender_fig1_20260923/PROVENANCE.json"
SURVEY = ROOT / "external/lit_survey_2026-09-13/lit_survey_background_controls.csv"
TALLY = ROOT / "results/claim_verdict_tally_2026-09-17.csv"
SKILL = pathlib.Path.home() / ".claude/skills/scipilot-figure-skill/scripts"
SELF = pathlib.Path(__file__).read_bytes()


def load(path, name):
    path = pathlib.Path(path)
    src = path.read_text(encoding="utf-8")
    if "__main__" not in src and re.search(r"^\s*(main|run)\(\)", src, re.M):
        sys.exit("REFUSE: %s runs on import" % path)
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m, hashlib.sha256(path.read_bytes()).hexdigest()


S, S_SHA = load(ROOT / "scripts/nc_figure_style_v2_2026-09-23.py", "style_v2")
CHK, CHK_SHA = load(ROOT / "scripts/check_figure_text_overlaps_2026-09-19.py", "overlap")
OLD, OLD_SHA = load(ROOT / "scripts/plot_figs_2to5_nc_style_2026-09-21.py", "old_gates")
VQA, VQA_SHA = load(SKILL / "visual_qa.py", "skill_vqa")
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402
import numpy as np  # noqa: E402

S.apply()
plt.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": "Arial"})

# step x-extent in the illustration, measured from its ink columns (pixel < 215), 2752 px wide -- unchanged from R35
STEPS = [  # (step, artefact number [colour key only, not drawn], artefact name, control as worded in R34/R35, ink span px)
    ("Enrichment", 4, "Positive set =\nidentified set",
     "Report the share of\nidentified cysteines with\na site; add a blocking arm", (58, 335)),
    ("Digestion", 1, "Cleavage\ngeometry", "Background matched\non theoretical\ndetectability", (533, 920)),
    ("LC–MS/MS", 3, "Protein\nabundance", "Abundance-matched\nprotein background", (1105, 1623)),
    ("Database search", 5, "Search\nspace", "Both candidate\nchemistries in the\nsearch space", (1798, 2187)),
    ("Site table", 2, "Multi-cysteine\npeptides", "Filter both sets to\nsingle-cysteine\npeptides", (2360, 2694)),
]
CROP = (20, 505, 2732, 1000)   # x0, y0, x1, y1 around the ink (rows 531-976, columns 58-2694) -- unchanged from R35
VLABEL = {"survives": "Survives", "undecidable": "Undecidable", "null_broken_by_control": "Null broken by control",
          "attenuated": "Attenuated", "vanishes": "Vanishes", "baseline_contradicts_claim": "Baseline contradicts claim"}


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def counts():
    rows = list(csv.DictReader(SURVEY.open(encoding="utf-8", newline="")))
    funnel = [("Records screened", sum(r["screen_included"] == "1" for r in rows)),
              ("Full text retrieved", sum(r["fulltext_retrieved"] == "1" for r in rows)),
              ("Sampled for coding", sum(r["in_classified_sample"] == "1" for r in rows)),
              ("Classified", sum(bool(r["background_class"].strip()) for r in rows))]
    cls = Counter(r["background_class"] for r in rows if r["background_class"].strip())
    comp = [("Other or not stated", cls["e"], S.GREY), ("All residues, unrestricted", cls["a"], S.GREY),
            ("Detected but unmodified", cls["c"], S.BLUE), ("Detectability or\nabundance matched", cls["b"], S.BLUE)]
    tally = [r for r in csv.DictReader(TALLY.open(encoding="utf-8", newline="")) if r["has_measured_baseline"] == "1"]
    verd = Counter(r["verdict"] for r in tally)
    order = ["survives", "undecidable", "null_broken_by_control", "attenuated", "vanishes", "baseline_contradicts_claim"]
    verdicts = [(VLABEL[v], verd[v], S.VERDICT[v]) for v in order]
    # G1: the values the manuscript already prints; a changed table refuses instead of redrawing
    assert [n for _, n in funnel] == [942, 367, 120, 74], funnel
    assert (cls["a"], cls["e"], cls["c"], cls["b"]) == (34, 35, 4, 1), cls
    assert [n for _, n, _ in verdicts] == [12, 11, 2, 1, 1, 1] and len(tally) == 28, verdicts
    return funnel, comp, verdicts, len(tally)


def hbars(ax, items, total_label=None, xmax=None):
    ys = list(range(len(items)))[::-1]
    for y, (lab, n, c) in zip(ys, items):
        ax.barh(y, n, height=0.62, color=c, lw=0, zorder=2)
        ax.text(n + (xmax or 1) * 0.02, y, str(n), va="center", ha="left", color=S.MUTE)
    ax.set_yticks(ys); ax.set_yticklabels([lab for lab, _, _ in items]); S.clean_y(ax)
    ax.set_ylim(-0.6, len(items) - 0.4)
    ax.set_xlim(0, xmax); ax.tick_params(axis="x", length=2)


def build():
    funnel, comp, verdicts, n_claims = counts()
    FIG_H = 110
    fig = plt.figure(figsize=(S.DOUBLE, FIG_H * S.MM))
    # --- a: illustration + annotation rows
    img = np.asarray(Image.open(ILLU).convert("RGB"))[CROP[1]:CROP[3], CROP[0]:CROP[2]].copy()
    img[(img >= 246).all(axis=2)] = 255
    h, w = img.shape[:2]
    left, right = 0.03, 0.985
    ax_img = fig.add_axes([left, 0.655, right - left, (right - left) * S.DOUBLE / (FIG_H * S.MM) * h / w])
    ax_img.imshow(img, interpolation="lanczos")
    ax_img.set_xlim(0, w); ax_img.set_ylim(h, 0); ax_img.axis("off")
    ax_img.set_xticks([]); ax_img.set_yticks([])
    top = ax_img.get_position().y1
    ax_note = fig.add_axes([left, 0.43, right - left, 0.21])
    ax_note.set_xlim(0, w); ax_note.set_ylim(0, 1); ax_note.axis("off")
    ax_note.set_xticks([]); ax_note.set_yticks([])
    for step, num, name, control, (x0, x1) in STEPS:
        cx = (x0 + x1) / 2 - CROP[0]
        fig.text(left + (right - left) * cx / w, top + 0.012, step, ha="center", va="bottom", fontsize=S.FS,
                 fontweight="bold", color=S.INK)
        colour = S.ARTEFACT[num]
        bx = x0 - CROP[0]  # left-aligned with the step's own ink, same anchor R35 used for the badge
        # R35 drew a coloured numbered badge here (FancyBboxPatch + digit) then the name in INK to its right.
        # R36 drops the badge and the digit; the artefact name itself carries the colour and sits at the badge's
        # old left edge, so the row keeps its identity through colour alone -- no number is drawn anywhere in panel a.
        ax_note.text(bx, 0.88, name, ha="left", va="center", color=colour, fontweight="bold",
                     fontsize=S.FS, linespacing=1.05)
        ax_note.text(bx, 0.60, control, ha="left", va="top", color=S.MUTE, fontsize=S.FS, linespacing=1.15)
    # --- b, c, d (unchanged from R35)
    gs = fig.add_gridspec(1, 3, left=0.17, right=0.975, bottom=0.105, top=0.36, wspace=1.05,
                          width_ratios=[1, 1, 1])
    ax_b = fig.add_subplot(gs[0])
    fcol = ["#c9d9e8", "#95b6d3", "#5f90bb", S.BLUE]
    hbars(ax_b, [(lab, n, c) for (lab, n), c in zip(funnel, fcol)], xmax=1150)
    ax_b.set_xticks([0, 500, 1000]); ax_b.set_xlabel("Records")
    ax_c = fig.add_subplot(gs[1])
    hbars(ax_c, comp, xmax=45)
    ax_c.set_xticks([0, 20, 40]); ax_c.set_xlabel("Classified analyses")
    ax_d = fig.add_subplot(gs[2])
    hbars(ax_d, verdicts, xmax=15)
    ax_d.set_xticks([0, 5, 10, 15]); ax_d.set_xlabel("Claims with a measured baseline")
    S.label_panels(fig, [(ax_b, "b", "Literature survey"), (ax_c, "c", "Background used"),
                         (ax_d, "d", "Re-tested claims (%d)" % n_claims)])
    fig.text(0.004, top + 0.012, "a", fontsize=S.FS_PANEL, fontweight="bold", ha="left", va="bottom", color=S.INK)
    return fig, dict(funnel=funnel, composition=[(l, n) for l, n, _ in comp],
                     verdicts=[(l, n) for l, n, _ in verdicts])


def main():
    final = "--draft" not in sys.argv
    prov = json.loads(PROV.read_text())
    if sha(ILLU) != prov["sha256"]:
        sys.exit("REFUSE: G6 illustration hash differs from PROVENANCE")
    outdir = ROOT / "figures_r36" if final else pathlib.Path(sys.argv[sys.argv.index("--draft") + 1])
    outdir.mkdir(parents=True, exist_ok=True)
    outs = [outdir / ("Fig1_overview.%s" % e) for e in ("pdf", "png", "svg")]
    aud = ROOT / "results/fig1_biorender_r36_2026-09-23_audit.json"
    if final and (any(p.exists() for p in outs) or aud.exists()):
        sys.exit("REFUSE: outputs exist")
    fig, values = build()
    ctl = plt.figure(); ctl.text(.5, .5, "ctl a", fontsize=10); ctl.text(.5, .5, "ctl b", fontsize=10)
    pc = CHK.audit_figure("ctl", ctl)["n_collisions"]; plt.close(ctl)
    g = OLD.gates_for(fig, "Fig1_overview")
    skill = VQA.audit_layout(fig)
    # G9: this round's whole point is that panel a no longer draws a digit. Read the rendered text objects in
    # ax_note (the annotation band, identified by its known axes position) and assert none of them is a bare 1-5.
    ax_note = [a for a in fig.axes if abs(a.get_position().y0 - 0.43) < 1e-6]
    assert len(ax_note) == 1, "could not locate ax_note for G9"
    digits_drawn = [t.get_text() for t in ax_note[0].texts if t.get_text().strip() in {"1", "2", "3", "4", "5"}]
    g9 = len(digits_drawn) == 0
    ok = (pc > 0 and g["collisions"] == 0 and g["min_lw"] >= 1.0 - 1e-9 and set(g["font_sizes"]) <= {7.0, 8.0}
          and not g["outside"] and not any(s == "FAIL" for s, _ in skill) and g9)
    print({k: g[k] for k in ("collisions", "pairs", "min_lw", "font_sizes", "outside", "size_mm")}, skill,
          {"g9_digits_drawn": digits_drawn})
    if ok:
        for p in outs:
            fig.savefig(p, dpi=600 if p.suffix == ".png" else None)
        VQA.render_preview(fig, str(outdir / "_preview_Fig1.png") if not final else
                           str(ROOT / "results/figure_skill_qa_r36_2026-09-23_Fig1_preview.png"), dpi=150)
    if final:
        aud.write_text(json.dumps(dict(
            script=pathlib.Path(__file__).name, script_sha256=hashlib.sha256(SELF).hexdigest(),
            imported=dict(style_v2=S_SHA, overlap=CHK_SHA, gates=OLD_SHA, skill_visual_qa=VQA_SHA),
            illustration=dict(path=str(ILLU.relative_to(ROOT)), sha256=sha(ILLU), crop=CROP,
                              provenance=str(PROV.relative_to(ROOT))),
            sources=dict(survey=str(SURVEY.relative_to(ROOT)), survey_sha256=sha(SURVEY),
                         tally=str(TALLY.relative_to(ROOT)), tally_sha256=sha(TALLY)),
            printed_values=values, gates=dict(positive_control=pc, g9_no_digits_in_panel_a=g9, **{k: g[k] for k in
                                              ("collisions", "min_lw", "font_sizes", "outside", "size_mm")}),
            skill_audit_layout=skill, outputs={str(p.relative_to(ROOT)): sha(p) for p in outs if p.exists()},
            note="panel a keeps identity through S.ARTEFACT colour only; no per-step number is drawn (user request "
                 "2026-09-23: renumbering the whole manuscript to workflow order was declined as too broad)",
            all_pass=ok), indent=2, ensure_ascii=False, default=str))
    print("all_pass", ok)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
