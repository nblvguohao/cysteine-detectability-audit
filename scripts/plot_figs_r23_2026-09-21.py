"""Build the R23 quantitative figure set.

Figs. 2–4 reuse the R22 source readers and stored values byte-for-byte. Fig. 5
keeps the same data and uncertainty definitions but changes the visual grammar:
panel b uses horizontal bars with asymmetric intervals and panel c uses a paired
profile across feature sets. This adds chart-type diversity without inventing
new estimands or subtracting rounded intervals.
"""
import csv
import hashlib
import importlib.util
import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "figures_r23"
SD = ROOT / "source_data_r23"
RES = ROOT / "results"
SELF = pathlib.Path(__file__).read_bytes()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BASE_PATH = ROOT / "scripts/plot_figs_2to5_nc_style_2026-09-21.py"
BASE = load(BASE_PATH, "r22_figure_base")
S = BASE.S
CHK = BASE.CHK
width_mm = 183


def read(path):
    return list(csv.DictReader(path.open(encoding="utf-8-sig", newline="")))


def panel_letter(ax, letter):
    ax.text(-0.02, 1.02, letter, transform=ax.transAxes, fontsize=S.FS_PANEL,
            fontweight="bold", ha="right", va="bottom")


def fig5():
    log = BASE.ValueLedger("Fig5")
    aud = read(RES / "self_audit_public_cohorts_2026-09-20.csv")
    abl = read(RES / "public_refit_ablation_2026-09-21.csv")
    five = read(RES / "self_audit_public_five_proteases_2026-09-20.csv")
    cohorts = [("human_PXD044043", "Human", BASE.HUMAN), ("rice_PXD072089", "Rice", BASE.RICE)]
    fig = plt.figure(figsize=(width_mm * S.MM, 128 * S.MM))
    gs = fig.add_gridspec(2, 2, left=0.205, right=0.985, bottom=0.105, top=0.95,
                          hspace=0.62, wspace=0.8)

    # a: retain the interval view because it is the primary evidence for
    # within-label-stratum ordering.
    ax = fig.add_subplot(gs[0, 0])
    y, ticks, labels = 0, [], []
    for coh, cn, col in cohorts:
        for st, sn in (("label_positive", "positives"), ("label_negative", "negatives")):
            for key, mn, c in (("v2_chem_deployed", "chemistry", col), ("v2_full_benchmark", "full", S.GREY)):
                row = next(x for x in aud if x["cohort"] == coh and x["model"] == key
                           and x["band"] == "distal_6_12" and x["stratum"] == st)
                S.ci(ax, y, float(row["auc"]), float(row["ci_low"]), float(row["ci_high"]), c)
                for field in ("auc", "ci_low", "ci_high"):
                    log.add("a", f"{coh}|{key}|{st}", field, row[field], float(row[field]),
                            "self_audit_public_cohorts_2026-09-20.csv")
                ticks.append(y); labels.append(f"{cn} {sn}, {mn}"); y -= 1
        y -= 0.6
    ax.axvline(0.5, color=S.GREY, lw=S.LW, ls=(0, (2, 2)))
    ax.set_yticks(ticks); ax.set_yticklabels(labels); ax.tick_params(axis="y", length=0)
    ax.set_ylim(y + 1.3, 0.7)
    ax.set_xlabel("AUC, score vs distal cleavage band")
    ax.set_title("Within each label stratum")
    panel_letter(ax, "a")

    # b: a bar/interval panel makes the control-versus-ablation comparison
    # visually distinct while preserving every CI in the source data.
    ax = fig.add_subplot(gs[0, 1])
    arms = [("control", "deployed"), ("drop_w710", "minus 2 windows"),
            ("drop_all_basic", "minus all basic"), ("label_permuted", "labels permuted")]
    y, ticks, labels = 0, [], []
    for coh, cn, col in cohorts:
        for arm, label in arms:
            row = next(x for x in abl if x["cohort"] == coh and x["arm"] == arm)
            point, lo, hi = (float(row[k]) for k in ("log2_or", "ci_low", "ci_high"))
            colour = S.GREY if arm == "label_permuted" else col
            ax.barh(y, point, height=0.52, color=colour, alpha=0.72,
                    edgecolor=colour, linewidth=S.LW, zorder=2)
            ax.errorbar(point, y, xerr=[[point - lo], [hi - point]], fmt="none",
                        ecolor=colour, elinewidth=S.LW, capsize=2, zorder=4)
            for field in ("log2_or", "ci_low", "ci_high"):
                log.add("b", f"{coh}|{arm}", field, row[field], float(row[field]),
                        "public_refit_ablation_2026-09-21.csv")
            ticks.append(y); labels.append(f"{cn}, {label}"); y -= 1
        y -= 0.6
    ax.axvline(0, color=S.INK, lw=S.LW)
    ax.set_yticks(ticks); ax.set_yticklabels(labels); ax.tick_params(axis="y", length=0)
    ax.set_ylim(y + 1.3, 0.7)
    ax.set_xlabel("Top-100 log2 odds ratio")
    ax.set_title("Deleting named columns")
    panel_letter(ax, "b")

    # c: paired profiles emphasize the pattern across arms rather than another
    # collection of horizontal dots.
    ax = fig.add_subplot(gs[1, 0])
    rows = [("control", "deployed"), ("drop_w710", "-2\nwindows"),
            ("drop_all_basic", "-all\nbasic"), ("detect_only", "detect.\nonly"),
            ("full", "all")]
    x = list(range(len(rows)))
    for coh, cn, col in cohorts:
        values = []
        for arm, _ in rows:
            row = next(r for r in abl if r["cohort"] == coh and r["arm"] == arm)
            values.append(float(row["within_protein_auc"]))
            log.add("c", f"{coh}|{arm}", "within_protein_auc", row["within_protein_auc"], values[-1],
                    "public_refit_ablation_2026-09-21.csv")
        ax.plot(x, values, marker="o", color=col, lw=S.LW, ms=4, label=cn, zorder=3)
    ax.set_xticks(x); ax.set_xticklabels([label for _, label in rows])
    ax.set_ylim(0.75, 0.92); ax.set_ylabel("Within-protein AUC")
    ax.set_title("Feature-set profile")
    ax.grid(axis="y", color="#E2E6E9", linewidth=S.LW, zorder=0)
    ax.legend(loc="lower left", ncol=2, frameon=False)
    panel_letter(ax, "c")

    # d: the cleavage-rule boundary remains an interval plot.
    ax = fig.add_subplot(gs[1, 1])
    proteases = ["Trypsin", "LysC", "GluC", "AspN", "Chymotrypsin"]
    for off, (coh, cn, col) in zip((-0.15, 0.15), cohorts):
        for i, protease in enumerate(proteases):
            row = next(r for r in five if r["cohort"] == coh and r["protease"] == protease
                       and r["model"] == "v2_chem_deployed" and r["band"] == "distal_6_12"
                       and r["stratum"] == "label_positive")
            S.ci(ax, i + off, float(row["auc"]), float(row["ci_low"]), float(row["ci_high"]), col)
            for field in ("auc", "ci_low", "ci_high"):
                log.add("d", f"{coh}|{protease}", field, row[field], float(row[field]),
                        "self_audit_public_five_proteases_2026-09-20.csv")
    ax.axvline(0.5, color=S.GREY, lw=S.LW, ls=(0, (2, 2)))
    ax.set_yticks(range(len(proteases))); ax.set_yticklabels(proteases); ax.invert_yaxis()
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("AUC, label-positive stratum")
    ax.set_title("Five cleavage rules")
    panel_letter(ax, "d")
    return fig, log


def render(name, result):
    fig, log = result[0], result[1]
    BASE.require_matplotlib_panel_alignment(fig)
    gate = BASE.gates_for(fig, name)
    OUT.mkdir(exist_ok=True); SD.mkdir(exist_ok=True)
    for suffix, kwargs in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600}), ("tiff", {"dpi": 600})):
        fig.savefig(OUT / f"{name}.{suffix}", **kwargs)
    plt.close(fig)
    with (SD / f"Source_Data_{name}.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(log.rows[0]))
        writer.writeheader(); writer.writerows(log.rows)
    gate["values_logged"] = len(log.rows)
    gate["pass"] = (gate["collisions"] == 0 and gate["min_lw"] >= 1.0
                     and set(gate["font_sizes"]).issubset({7.0, 8.0}) and not gate["outside"])
    return gate


def main():
    OUT.mkdir(exist_ok=True); SD.mkdir(exist_ok=True)
    control = plt.figure()
    control.text(.5, .5, "ctl", fontsize=10); control.text(.5, .5, "ctl", fontsize=10)
    positive_control = CHK.audit_figure("positive_control", control)["n_collisions"]
    plt.close(control)
    figures = {}
    # Re-render Figs. 2–4 using the same stored readers, but into R23 folders.
    for name, fn in (("Fig2_search_space", BASE.fig2), ("Fig3_three_axes", BASE.fig3),
                     ("Fig4_claim_retests", BASE.fig4), ("Fig5_self_audit_public", fig5)):
        figures[name] = render(name, fn())
        print(name, json.dumps(figures[name], sort_keys=True))
    ok = positive_control > 0 and all(x["pass"] for x in figures.values())
    audit = {
        "script": str(pathlib.Path(__file__).relative_to(ROOT)),
        "script_sha256": hashlib.sha256(SELF).hexdigest(),
        "base_script_sha256": hashlib.sha256(BASE_PATH.read_bytes()).hexdigest(),
        "positive_control": positive_control,
        "figures": figures,
        "fig5_design_change": {
            "panel_b": "horizontal bars with asymmetric intervals",
            "panel_c": "paired feature-set profiles",
            "data_change": False,
            "uncertainty_change": False,
        },
        "all_pass": ok,
    }
    (RES / "figs_r23_2026-09-21_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
