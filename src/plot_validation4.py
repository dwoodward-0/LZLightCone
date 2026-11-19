"""Validation 4 figure: total background (all 12 sources) vs WS2024 Fig. 4.

Draws the declared background_projection_validation figure. It runs AFTER
src/validate_backgrounds.py in the same recipe. That script evaluates the six
acceptance gates and writes them to a JSON sidecar next to --output (and a
per-category diagnostic PNG at --output, which this script replaces). We read
the gates from that JSON and print them as the figure footer, so the figure
and the machine-readable verdict always come from the same run.

The figure itself is src/plot_background_total.draw_total (formerly the
review figure results/bkg_review/bkg_total.png; made THE Validation 4 figure
by user decision 2026-11-19):
  - black step: LZ Fig. 4 drawn best-fit total background;
  - blue band: Fig. 4's displayed uncertainty envelope (systematic + stat.),
    digitized directly from the pinned PDF's vector fills;
  - red dashed step and band: our sum of all 12 Table I sources, +-1 sigma
    finite-MC statistics only.
Presentation convention (user decision 2026-11-11): the red curve shows 124Xe
at the Fig.-4 shape-optimal LL charge yield q = 0.78 (analytic S2 rescaling of
S1c > 40 events; LZ's fitted value is 0.70 +- 0.04). The stored model and the
six gates use the production q = 0.65, except gate 5's 124Xe entry, which is
evaluated at the best q within LZ's prior. The footer says so.
"""
import argparse
import json
from pathlib import Path

from plot_background_total import draw_total
from plot_backgrounds_per_source import XE124_Q_LL_BEST_RMS

# Gate thresholds, mirrored from src/validate_backgrounds.py (its docstring
# and main()); used only to label the footer.
PANEL_TOL = 0.20
TOTAL_RMS_MAX = 0.15
CATEGORY_RMS_MAX = 0.25
TAIL_BUDGET = 8.0


def footer_text(g):
    acc = g["acceptance"]
    n_pass = sum(bool(acc[k]) for k in
                 ("total", "components", "panels", "total_shape",
                  "categories", "tail"))
    ok = lambda k: "PASS" if acc[k] else "FAIL"
    pt = g["panel_totals"]
    rel = {p: pt[p]["ours"] / pt[p]["lz_digitized"] - 1.0 for p in pt}
    xe = g["xe124_charge_yield"]
    cats = g["residuals_dex"]
    cat_rms = max(v["rms"] for k, v in cats.items()
                  if k in ("pb214", "other_er_as_drawn", "xe127",
                           "accidentals") and v)
    cat_rms = max(cat_rms, xe["rms_at_best"])
    lines = [
        f"Validation 4: {'PASS' if acc['pass'] else 'FAIL'} "
        f"({n_pass}/6 gates)",
        f"1 total {g['total_counts']:.1f} vs Table I {g['table1_total']:.0f}"
        f" +/- 91 [{ok('total')}]   2 per-source = Table I [{ok('components')}]",
        f"3 panel totals top {rel['top']:+.1%}, bottom {rel['bottom']:+.1%}"
        f" (<= {PANEL_TOL:.0%}) [{ok('panels')}]   4 total-curve rms "
        f"{cats['total']['rms']:.3f} dex (<= {TOTAL_RMS_MAX}) "
        f"[{ok('total_shape')}]",
        f"5 per-category rms max {cat_rms:.3f} dex (<= {CATEGORY_RMS_MAX})"
        f" [{ok('categories')}]   6 sub-floor tail "
        f"{g['tail_abs_events']['total']:.1f} events (<= {TAIL_BUDGET:.0f})"
        f" [{ok('tail')}]",
        f"Gates use the stored model (124Xe LL q = {xe['pre_fit']:.2f}; gate 5"
        f" 124Xe at best q = {xe['best_rms_q']:.2f}).",
        f"Red curve: 124Xe at Fig.-4 shape-optimal q = "
        f"{XE124_Q_LL_BEST_RMS:.2f} (LZ fit 0.70 +/- 0.04).",
    ]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--models", required=True)
    ap.add_argument("--response", required=True)
    ap.add_argument("--paper", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    gates_path = Path(args.output).with_suffix(".json")
    gates = json.loads(gates_path.read_text())
    footer = footer_text(gates)
    print(footer)
    sums = draw_total(args.models, args.response, args.paper, args.output,
                      footer=footer)
    gates["figure_total_curve_sums_q078"] = sums
    gates_path.write_text(json.dumps(gates, indent=2))


if __name__ == "__main__":
    main()
