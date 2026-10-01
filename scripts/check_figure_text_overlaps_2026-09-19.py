"""Programmatic text-collision audit for the four submission figures.

WHY THIS EXISTS
The user reported seeing text-on-text and text-on-table-row overlaps in earlier renders of
this project's figures. I had only checked collisions by eye on downsampled screenshots
(CLAUDE.md 9.56.3 items 6-9 record four defects found that way). Eyeballing a screenshot is
exactly the failure mode this tree's house rules exist to prevent: it depends on me noticing,
and a screenshot rescaled for display can hide a few-pixel overlap that is real in the 600 dpi
file. This script measures instead.

WHAT IT DOES
For each of F1-F4, in their final journalised state (bold panel letters, print-resolution
sizing), it collects every Text artist actually drawn: axis titles, x/y labels, tick labels,
legend entries and legend titles, and every free-floating ax.text/annotate call. It renders the
canvas once so extents are real (not the pre-layout estimate), converts every text's tight
bounding box to display (pixel) coordinates, and checks all pairs for a strictly positive
intersection area. There are no matplotlib Table artists in these figures (grep confirms it);
row-label lists and value columns are drawn as ordinary Text, so they are already covered.

A small, explicit tolerance (2 px in the print-resolution frame) is subtracted from each box
before intersection, because adjacent tick labels legitimately share a boundary pixel; without
it every dense axis would flag as "colliding" and the report would be useless. Anything larger
than that tolerance is reported by the two texts' content, not just coordinates, so a human can
tell instantly what to fix.

WHY THIS SCRIPT DOES NOT EDIT THE PLOT SCRIPTS OR THE RENDER WRAPPER
scripts/plot_manuscript_f{1,2,3,4}*.py are referenced by sha256 elsewhere and are byte-identical
per CLAUDE.md 9.56.4. scripts/render_figures_submission_2026-09-19.py is likewise already cited
by sha256 in results/figures_submission_render_2026-09-19_audit.json and in CLAUDE.md 9.56. Both
are import-safe (guarded by __main__), so this script imports the render module and MONKEY-PATCHES
ITS OWN `patched_savefig` NAME to intercept the fully journalised Figure object right before it is
written and closed, then calls its unmodified main(). Neither file is edited.

INTERPRETER: /Users/lyuguohao/Documents/lnrna复现/.venv/bin/python (matplotlib 3.9.4), the one
CLAUDE.md 9.20.7 fixes for every figure in this tree.
"""
import hashlib, importlib.util, json, pathlib, sys
import matplotlib
matplotlib.use("Agg")

ROOT = pathlib.Path(__file__).resolve().parent.parent
SELF = pathlib.Path(__file__).read_bytes()
RENDER_SCRIPT = ROOT / "scripts/render_figures_submission_2026-09-19.py"

TOL_PX = 2.0   # adjacent labels may share this many pixels of border without being a defect

def collect_texts(fig):
    """Every Text artist actually on the canvas: titles (all three loc slots), axis labels,
    tick labels, legend entries/titles, and free ax.text/annotate calls. Empty strings and
    invisible artists are skipped."""
    out = []
    for ax in fig.get_axes():
        cands = []
        for loc_attr in ("_left_title", "title", "_right_title"):
            t = getattr(ax, loc_attr, None)
            if t is not None:
                cands.append(t)
        cands += [ax.xaxis.label, ax.yaxis.label]
        cands += list(ax.get_xticklabels()) + list(ax.get_yticklabels())
        cands += list(ax.texts)
        leg = ax.get_legend()
        if leg is not None:
            cands += list(leg.get_texts())
            if leg.get_title() is not None:
                cands.append(leg.get_title())
        for t in cands:
            if t.get_text().strip() and t.get_visible():
                out.append(t)
    for t in fig.texts:
        if t.get_text().strip() and t.get_visible():
            out.append(t)
    return out


def bbox_of(text, renderer):
    b = text.get_window_extent(renderer=renderer)
    return (b.x0 + TOL_PX, b.y0 + TOL_PX, b.x1 - TOL_PX, b.y1 - TOL_PX)


def overlaps(a, b):
    ax0, ay0, ax1, ay1 = a
    bx0, by0, bx1, by1 = b
    ix0, iy0 = max(ax0, bx0), max(ay0, by0)
    ix1, iy1 = min(ax1, bx1), min(ay1, by1)
    return max(0.0, ix1 - ix0) * max(0.0, iy1 - iy0)


def audit_figure(stem, fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = collect_texts(fig)
    boxes = [bbox_of(t, renderer) for t in texts]
    findings = []
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            area = overlaps(boxes[i], boxes[j])
            if area > 0:
                findings.append(dict(
                    text_a=texts[i].get_text()[:60], text_b=texts[j].get_text()[:60],
                    overlap_area_px2=round(area, 1),
                    a_bbox=[round(v, 1) for v in boxes[i]], b_bbox=[round(v, 1) for v in boxes[j]]))
    return dict(figure=stem, n_texts=len(texts), n_collisions=len(findings), collisions=findings)


def main():
    src = RENDER_SCRIPT.read_text()
    assert '__main__' in src, "render script is not import-safe; refusing to import"
    render_sha = hashlib.sha256(RENDER_SCRIPT.read_bytes()).hexdigest()

    spec = importlib.util.spec_from_file_location("render_mod", RENDER_SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)   # safe: guarded, does not auto-run main()

    captured = []
    original_patched_savefig = m.patched_savefig

    def capturing_savefig(self, fname, **kw):
        if str(fname).endswith(".png"):
            stem = pathlib.Path(str(fname)).stem
            captured.append((stem, self))
        return original_patched_savefig(self, fname, **kw)

    m.patched_savefig = capturing_savefig   # m.main() reads this name from its own globals
    m.main()

    results = [audit_figure(stem, fig) for stem, fig in captured]
    total_collisions = sum(r["n_collisions"] for r in results)
    gates = {
        "G1_all_four_figures_captured": sorted(r["figure"].split("_")[1] for r in results) == ["f1", "f2", "f3", "f4"],
        "G2_zero_text_collisions": total_collisions == 0,
        "G2_total_collisions": total_collisions,
        "G3_lookup_found_something": sum(r["n_texts"] for r in results) > 0,
    }
    (ROOT/"results"/"figure_text_overlap_audit_2026-09-19.json").write_text(json.dumps(dict(
        script=pathlib.Path(__file__).name, script_sha256=hashlib.sha256(SELF).hexdigest(),
        imported_render_script=dict(name=RENDER_SCRIPT.name, sha256=render_sha, import_safe_guard_present=True),
        tolerance_px=TOL_PX, per_figure=results, gates=gates,
        what_this_does_not_check=[
            "Text-vs-line or text-vs-marker collisions (error bars, dashed reference lines) are "
            "not checked; the user's report named text-vs-text and text-vs-table specifically.",
            "It checks the four submission-spec renders, not the un-journalised originals in results/."]),
        indent=2, ensure_ascii=False))
    for r in results:
        print(f"  {r['figure']:44s} texts={r['n_texts']:3d}  collisions={r['n_collisions']}")
        for c in r["collisions"]:
            print(f"      COLLIDE  {c['overlap_area_px2']:7.1f}px2  {c['text_a']!r}  <->  {c['text_b']!r}")
    print("\nGATES:", {k: v for k, v in gates.items() if "collisions" not in results})
    if not (gates["G1_all_four_figures_captured"] and gates["G3_lookup_found_something"]):
        sys.exit(1)
    print("\n" + ("ALL CLEAR: 0 text collisions across all four figures" if gates["G2_zero_text_collisions"]
          else f"FOUND {total_collisions} collision(s) - see above, not auto-fixed by this script"))

if __name__ == "__main__":
    main()
