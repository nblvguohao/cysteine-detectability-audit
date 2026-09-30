"""Render the R23 overview figure for the Nature Communications manuscript.

Figure contract
---------------
Core conclusion: site-level modification claims become interpretable only when
the control matches the workflow property that generated each observable set.
Results-level question: where do the five artefacts enter a bottom-up workflow,
and which control makes the corresponding inference defensible?
Archetype: schematic-led composite, with a hero measurement-to-inference map
and a compact control/readout matrix.
Backend: Python/Matplotlib, drawn at the 183-mm double-column width.

The figure is conceptual. It introduces no numerical result; quantitative
results remain in Figs. 2–5 and their source-data files. Labels are deliberately
short so the figure carries the visual argument while the legend carries the
definitions. The five colours are the manuscript's fixed artefact mapping.
"""
import hashlib
import importlib.util
import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle, FancyArrowPatch


ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "figures_r23"
SELF = pathlib.Path(__file__).read_bytes()


def load(name):
    path = ROOT / "scripts" / name
    spec = importlib.util.spec_from_file_location(name[:-3].replace("-", "_"), path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, hashlib.sha256(path.read_bytes()).hexdigest()


STYLE, STYLE_SHA = load("nc_figure_style_2026-09-21.py")
CHECK, CHECK_SHA = load("check_figure_text_overlaps_2026-09-19.py")
STYLE.apply()
matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})

NAVY = "#17324D"
PALE = "#F5F7F8"
PALE_BLUE = "#EAF1F6"
MID = "#AAB5BE"
WHITE = "#FFFFFF"

STAGES = [
    ("Enrichment", "capture / release"),
    ("Digestion", "protease"),
    ("LC–MS/MS", "selection"),
    ("Search", "candidate space"),
    ("Site table", "assigned residues"),
    ("Claim test", "comparison"),
]

ARTEFACTS = [
    (1, "Positive ≈ identified", "no observed\nnegative class", 0),
    (2, "Cleavage geometry", "K/R spacing\nmimics motif", 1),
    (3, "Protein abundance", "detection tracks\nabundance", 2),
    (4, "Search space", "omitted chemistry\ncannot be called", 3),
    (5, "Multi-Cys peptides", "site assignment\ncan flip signs", 4),
]

MATRIX = [
    (1, "Observed negative arm", "negative class exists"),
    (2, "Detectability-matched background", "sequence signal is testable"),
    (3, "Abundance-matched proteins", "protein commonality is interpretable"),
    (4, "Both candidate chemistries", "chemical identity is resolvable"),
    (5, "Single-Cys filtering in both sets", "localisation sign is stable"),
]

OUTPUTS = [
    ("CHECKLIST", "pre-analysis"),
    ("LABEL TIERS", "T1 → T3"),
    ("BENCHMARK", "scope + limits"),
    ("SELF-AUDIT", "model on public data"),
]


def rounded(ax, x, y, w, h, fc=WHITE, ec=STYLE.INK, lw=STYLE.LW, radius=1.5, z=1):
    patch = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={radius}",
                           facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z)
    ax.add_patch(patch)
    return patch


def arrow(ax, x0, y0, x1, y1, color=STYLE.INK, lw=STYLE.LW, ms=7, z=5):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=ms, linewidth=lw, color=color,
                                 shrinkA=0, shrinkB=0, zorder=z))


def badge(ax, x, y, number, radius=2.6, z=7):
    ax.add_patch(Circle((x, y), radius, facecolor=STYLE.ARTEFACT[number],
                        edgecolor=WHITE, linewidth=STYLE.LW, zorder=z))
    ax.text(x, y, str(number), ha="center", va="center", color=WHITE,
            fontsize=STYLE.FS, fontweight="bold", zorder=z + 1)


def text(ax, x, y, value, **kwargs):
    defaults = dict(fontsize=STYLE.FS, color=STYLE.INK, zorder=8)
    defaults.update(kwargs)
    return ax.text(x, y, value, **defaults)


def main():
    width_mm, height_mm = 183, 126
    fig = plt.figure(figsize=(width_mm * STYLE.MM, height_mm * STYLE.MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width_mm)
    ax.set_ylim(0, height_mm)
    ax.axis("off")

    # Panel a: a measurement-to-inference map. The assay is a neutral spine;
    # coloured risk pins mark where an artefact enters the comparison.
    text(ax, 2, 124, "a", fontsize=STYLE.FS_PANEL, fontweight="bold", va="top")
    text(ax, 8, 124, "From an observed cysteine to a site-level claim", va="top", fontweight="bold")
    text(ax, 8, 119, "The assay produces an observable set; the audit asks whether the comparison preserves its meaning.",
         color=STYLE.GREY)

    left, y_stage, stage_w, gap, stage_h = 7, 103, 24, 4, 13
    stage_centres = []
    for i, (title, subtitle) in enumerate(STAGES):
        x = left + i * (stage_w + gap)
        stage_centres.append(x + stage_w / 2)
        final = i == len(STAGES) - 1
        rounded(ax, x, y_stage, stage_w, stage_h, fc=NAVY if final else PALE,
                ec=NAVY if final else MID, lw=STYLE.LW, radius=1.8)
        text(ax, x + stage_w / 2, y_stage + 8.1, title, ha="center", va="center",
             color=WHITE if final else STYLE.INK, fontweight="bold")
        text(ax, x + stage_w / 2, y_stage + 3.5, subtitle, ha="center", va="center",
             color=WHITE if final else STYLE.GREY)
        if i < len(STAGES) - 1:
            arrow(ax, x + stage_w + 0.8, y_stage + stage_h / 2,
                  x + stage_w + gap - 0.8, y_stage + stage_h / 2, color=MID, ms=5)

    # Five compact callouts. The small pins establish a visual link to the
    # assay stage without creating a second long flowchart.
    card_y, card_h, card_w = 78, 17, 32
    card_xs = [7, 43, 79, 115, 151]
    for (number, title, detail, stage_index), x in zip(ARTEFACTS, card_xs):
        target_x = stage_centres[stage_index]
        badge(ax, target_x, y_stage - 3.3, number, radius=2.4)
        ax.plot([target_x, target_x], [y_stage - 5.7, card_y + card_h],
                color=STYLE.ARTEFACT[number], linewidth=STYLE.LW, zorder=2)
        rounded(ax, x, card_y, card_w, card_h, fc=WHITE, ec=STYLE.ARTEFACT[number],
                lw=STYLE.LW, radius=1.8)
        ax.add_patch(Rectangle((x, card_y), 2.6, card_h, facecolor=STYLE.ARTEFACT[number],
                               edgecolor="none", zorder=2))
        text(ax, x + 5, card_y + 11.5, title, va="center", fontweight="bold")
        text(ax, x + 5, card_y + 5.2, detail, va="center", linespacing=1.15)

    # The inferential gate is the hero element of panel a.
    rounded(ax, 7, 58, 169, 13, fc=PALE_BLUE, ec=NAVY, lw=STYLE.LW, radius=2.0)
    text(ax, 12, 64.5, "OBSERVABLE SET", color=NAVY, fontweight="bold", va="center")
    arrow(ax, 47, 64.5, 61, 64.5, color=NAVY, ms=6)
    text(ax, 66, 64.5, "one control matched to one artefact", color=NAVY, fontweight="bold", va="center")
    arrow(ax, 123, 64.5, 136, 64.5, color=NAVY, ms=6)
    rounded(ax, 140, 60.2, 32, 8.6, fc=NAVY, ec=NAVY, lw=STYLE.LW, radius=1.5)
    text(ax, 156, 64.5, "CLAIM", color=WHITE, fontweight="bold", ha="center", va="center")
    text(ax, 91.5, 54.5, "If the control is absent, the result may describe the workflow rather than the chemistry.",
         ha="center", va="top", color=STYLE.GREY, style="italic")

    # Panel b: matrix plus four reusable outputs. This is intentionally not a
    # second flowchart: it is a compact reading aid for the five controls.
    text(ax, 2, 49, "b", fontsize=STYLE.FS_PANEL, fontweight="bold", va="top")
    text(ax, 8, 49, "One control per artefact, four reusable outputs", va="top", fontweight="bold")
    text(ax, 8, 43.5, "AUDIT QUESTION", color=STYLE.GREY, fontweight="bold")
    text(ax, 71, 43.5, "CONTROL", color=STYLE.GREY, fontweight="bold")
    text(ax, 137, 43.5, "READOUT", color=STYLE.GREY, fontweight="bold")
    ax.plot([7, 176], [42, 42], color=MID, linewidth=STYLE.LW)

    row_y = [37.5, 31, 24.5, 18, 11.5]
    for (number, control, readout), y in zip(MATRIX, row_y):
        ax.add_patch(Rectangle((7, y - 2.3), 2.2, 4.6, facecolor=STYLE.ARTEFACT[number],
                               edgecolor="none", zorder=2))
        badge(ax, 13, y, number, radius=2.0)
        text(ax, 18, y, ARTEFACTS[number - 1][1], va="center", fontweight="bold")
        text(ax, 71, y, control, va="center")
        text(ax, 137, y, readout, va="center")
        if y != row_y[-1]:
            ax.plot([7, 176], [y - 3.2, y - 3.2], color="#E2E6E9", linewidth=STYLE.LW)

    output_y, output_h, output_w, output_gap = 0.8, 7.0, 33, 2.2
    text(ax, 7, 4.3, "REUSABLE OUTPUTS", color=STYLE.GREY, fontweight="bold")
    for i, (title, subtitle) in enumerate(OUTPUTS):
        x = 39 + i * (output_w + output_gap)
        rounded(ax, x, output_y, output_w, output_h, fc=NAVY if i == 0 else PALE,
                ec=NAVY if i == 0 else MID, lw=STYLE.LW, radius=1.4)
        text(ax, x + 3.5, output_y + 4.0, title, va="center", fontweight="bold",
             color=WHITE if i == 0 else NAVY)
        text(ax, x + 3.5, output_y + 1.8, subtitle, va="center",
             color=WHITE if i == 0 else STYLE.GREY)

    # Render-time QA. Number gate is intentionally empty: the schematic does
    # not introduce quantitative values.
    fig.canvas.draw()
    render = fig.canvas.get_renderer()
    fw, fh = fig.canvas.get_width_height()
    outside = []
    for artist in [*ax.texts, *fig.texts]:
        if not (artist.get_text() or "").strip():
            continue
        box = artist.get_window_extent(render)
        if box.x0 < -0.5 or box.y0 < -0.5 or box.x1 > fw + 0.5 or box.y1 > fh + 0.5:
            outside.append(artist.get_text()[:60])
    collision = CHECK.audit_figure("fig1_overview_r23", fig)
    control = plt.figure()
    control.text(0.5, 0.5, "ctl", fontsize=10)
    control.text(0.5, 0.5, "ctl", fontsize=10)
    positive_control = CHECK.audit_figure("positive_control", control)["n_collisions"]
    plt.close(control)
    gates = {
        "G1_no_new_numeric_claims": True,
        "G2_min_linewidth": STYLE.min_linewidth(fig),
        "G3_font_sizes": STYLE.font_sizes(fig),
        "G4_text_collisions": collision["n_collisions"],
        "G4_positive_control": positive_control,
        "G5_text_outside_canvas": outside,
        "G6_dimensions_mm": [width_mm, height_mm],
    }
    OUT.mkdir(exist_ok=True)
    for suffix, kwargs in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600}), ("tiff", {"dpi": 600})):
        fig.savefig(OUT / f"Fig1_overview.{suffix}", **kwargs)
    plt.close(fig)
    ok = (gates["G2_min_linewidth"] >= 1.0 and set(gates["G3_font_sizes"]).issubset({7.0, 8.0})
          and gates["G4_text_collisions"] == 0 and gates["G4_positive_control"] > 0
          and not gates["G5_text_outside_canvas"])
    audit = {
        "script": str(pathlib.Path(__file__).relative_to(ROOT)),
        "script_sha256": hashlib.sha256(SELF).hexdigest(),
        "style_sha256": STYLE_SHA,
        "detector_sha256": CHECK_SHA,
        "size_mm": [width_mm, height_mm],
        "archetype": "schematic-led composite",
        "figure_contract": {
            "core_conclusion": "site-level claims require artefact-matched controls",
            "panel_a": "measurement-to-inference map",
            "panel_b": "control/readout matrix and reusable outputs",
        },
        "gates": gates,
        "all_pass": ok,
    }
    (ROOT / "results/fig1_overview_r23_2026-09-21_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
