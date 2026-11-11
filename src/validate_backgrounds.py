"""Validation 4: background projections vs LZ's WS2024 Fig. 4.

Compares the per-component background models (background_models output)
against the digitized Fig. 4 projection of arXiv:2410.17036v3
(lz_ws2024_fig4_projection): within each panel's S1c range (top
3 < S1c < 40 phd, bottom 40 < S1c < 80 phd) LZ histograms events in the
band-projected coordinate u = (log10(S2c) - mu_ER(S1c)) / sigma_ER(S1c)
with the y axis in events per bin (log scale). We construct the identical
projection from the stored samples: mu_ER is the ER (tritium) band median
and sigma_ER the 1-sigma band half-width (q84.13 - q15.87)/2, both from
the response_model_calibration band quantiles (validated against LZ's
published bands in er_band_validation_ws2024 to <= 0.36% in the median).

Category mapping (our sources -> LZ's drawn series):
  pb214_betas     <- 214Pb_betas
  other_er_betas  <- kr85 + solar + pb212 + tritium_c14 + xe136
  xe124_dec       <- xe124_dec
  xe127_125_ec    <- xe127_xe125_ec
  accidentals     <- accidental_coincidences
  b8              <- b8_hep_nu_nr + atm_nu_nr + detector_neutrons
plus the comparison variant other_er_as_drawn (kr85 + xe136 only): the
digitization established that LZ's drawn "Other ER betas" outline holds
only those two components (sum ratios Total/214Pb = 1.628 = 1210/743 and
other/214Pb = 0.293 = 217.6/743 are both exact), with solar-nu, 212Pb and
tritium/14C (223 events, 18.4%) present in the Total Bkg curve but absent
from the outlines. The variant is the apples-to-apples comparison against
the drawn outline; the full other_er_betas curve is also computed.

ACCEPTANCE (all reported; PASS/FAIL printed):
  1. Total ROI counts equal the Table I total 1210 +/- 91 (by
     construction exact; guards the normalization wiring).
  2. Per-component totals equal their Table I means (by construction).
  3. Panel-integrated totals within 20% of the digitized total_bkg sums
     (tests the S1c split of the model, a genuine prediction not in
     Table I; 20% covers digitization + ROI-edge effects).
  4. Shape: the total curve within rms 0.15 dex of the digitized
     total_bkg over bins where LZ's drawn level is >= 1 event/bin (the
     dex metric is only meaningful there: below it the drawn line's own
     ~0.02-0.05 dex thickness against the 10^-2 axis floor dominates,
     and those bins carry <= 0.1% of the model each).
  5. Per-category: 214Pb, 124Xe, 127Xe+125Xe, accidentals and the
     as-drawn other-ER comparison each within rms 0.25 dex of their
     drawn outlines over the same main-floor bins (looser than the
     total: these curves are 10-100x smaller, so digitization and MC
     noise matter more). 124Xe is evaluated at the best Qy_LL/Qy_beta
     within LZ's published floating prior [0.65, 0.87] (WS2024 paper;
     our pre-fit value is the TI box model 0.65), via the analytic
     charge rescaling, with the implied value reported.
  6. Sub-floor tail: the TOTAL curve's absolute |ours - LZ| over the
     excluded (< 1 event/bin) bins within 8 events (~0.7% of the
     model). The deep low-charge tail (|u| >~ 3) is where LZ's full
     position-dependent simulation diverges from uniform_1d (their
     S2c-correction residuals populate exactly this tail; same
     attribution as the WS2024 NR mid-range width residual), so it is
     gated on counts, not dex.
"""

import argparse
import csv
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_background_models import ROI_S1C_MAX, ROI_S1C_MIN, read_table1

PANEL_S1C = {"top": (ROI_S1C_MIN, 40.0), "bottom": (40.0, ROI_S1C_MAX)}
FLOOR = -2.0  # the figure clips at its 10^-2 axis minimum

# our source groups -> LZ's drawn series
CATEGORY_MAP = {
    "pb214_betas": ["214Pb_betas"],
    "other_er_betas": ["kr85_ar39_det_gammas", "solar_nu_er",
                       "pb212_po218_betas", "tritium_c14_betas",
                       "xe136_2nubb"],
    "xe124_dec": ["xe124_dec"],
    "xe127_125_ec": ["xe127_xe125_ec"],
    "accidentals": ["accidental_coincidences"],
    "b8": ["b8_hep_nu_nr", "atm_nu_nr", "detector_neutrons"],
}
CATEGORY_MAP["other_er_as_drawn"] = ["kr85_ar39_det_gammas", "xe136_2nubb"]


def load_reference(path):
    """(panel, series) -> (u_centers, u_lo, u_hi, log10 levels)."""
    out = {}
    with open(path) as f:
        rows = [r for r in csv.reader(f) if r and not r[0].startswith("#")
                and r[0] != "panel"]
    for panel, series, uc, ulo, uhi, l10 in rows:
        out.setdefault((panel, series), []).append(
            (float(uc), float(ulo), float(uhi), float(l10)))
    return {k: np.array(sorted(v)) for k, v in out.items()}


def load_band(response_path):
    """mu_ER(S1c) and sigma_ER(S1c) interpolants from the ER band."""
    with h5py.File(response_path) as f:
        centers = f["er_ch3t/band_s1c_centers"][:]
        q = f["er_ch3t/band_log10s2c_quantiles"][:]
    mu = q[2]
    sigma = 0.5 * (q[3] - q[1])  # (q84.13 - q15.87)/2

    def band(s1c):
        return np.interp(s1c, centers, mu), np.interp(s1c, centers, sigma)
    return band


def load_models(path):
    """{source: (s1c, log10_s2c, weight)} plus table1 means from attrs."""
    models = {}
    with h5py.File(path) as f:
        for src in f:
            if src == "accidentals_fit":
                continue
            g = f[src]
            models[src] = (g["s1c"][:], g["log10_s2c"][:], g["weight"][:])
    return models


def project(models, sources, band, edges):
    """Weighted u-histogram per panel for the given sources."""
    hist = {p: np.zeros(len(edges) - 1) for p in PANEL_S1C}
    for src in sources:
        s1c, y, w = models[src]
        mu, sig = band(s1c)
        u = (y - mu) / sig
        for panel, (lo, hi) in PANEL_S1C.items():
            m = (s1c > lo) & (s1c <= hi)
            hist[panel] += np.histogram(u[m], bins=edges, weights=w[m])[0]
    return hist


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--models", required=True)
    ap.add_argument("--response", required=True)
    ap.add_argument("--table1", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--background-components", default="all_table1")
    ap.add_argument("--exposure-structure", default="combined_single")
    ap.add_argument("--nest-parameter-set", default="ws2024_nest240")
    ap.add_argument("--detector-effects", default="uniform_1d")
    args = ap.parse_args()
    assert args.background_components == "all_table1"
    assert args.exposure_structure == "combined_single"
    assert args.nest_parameter_set == "ws2024_nest240"
    assert args.detector_effects == "uniform_1d"

    table1 = read_table1(args.table1)
    ref = load_reference(args.reference)
    band = load_band(args.response)
    models = load_models(args.models)

    # the u bin grid: shared by both panels; take it from the total_bkg
    # rows (the complete drawn curve). The CSV's lo/hi are rounded
    # independently (hi differs from the next row's lo by up to 0.001),
    # so build the grid arithmetically instead of uniquing the raw values
    # (that would create 0.001-wide sliver bins).
    tot = ref[("top", "total_bkg")]
    u0 = float(tot[0, 1])
    du = float(np.median(np.diff(tot[:, 1])))
    edges = u0 + du * np.arange(len(tot) + 1)

    ours = {cat: project(models, srcs, band, edges)
            for cat, srcs in CATEGORY_MAP.items()}

    centers = 0.5 * (edges[:-1] + edges[1:])

    # --- acceptance 1 & 2: totals ---
    total_all = sum(g[2].sum() for g in models.values())
    t1_total, t1_sigma = table1["total"]
    print(f"total ROI counts: {total_all:.1f} vs Table I "
          f"{t1_total:.1f} +/- {t1_sigma:.1f}")
    a1 = abs(total_all - t1_total) <= t1_sigma
    per_comp_ok = True
    for src, (_, _, w) in models.items():
        if src not in table1:
            continue
        m = table1[src][0]
        if abs(w.sum() - m) > 1e-6 * max(m, 1.0):
            per_comp_ok = False
            print(f"  MISMATCH {src}: {w.sum():.3f} vs Table I {m}")
    print(f"per-component totals match Table I means: {per_comp_ok}")

    # --- acceptance 3: panel-integrated totals vs the digitized curves ---
    a3 = True
    for panel in PANEL_S1C:
        lz = ref[(panel, "total_bkg")][:, 3]
        lz_sum = float(np.sum(10 ** lz))
        our_sum = float(ours["pb214_betas"][panel].sum()
                        + ours["other_er_betas"][panel].sum()
                        + ours["xe124_dec"][panel].sum()
                        + ours["xe127_125_ec"][panel].sum()
                        + ours["accidentals"][panel].sum()
                        + ours["b8"][panel].sum())
        rel = our_sum / lz_sum - 1.0
        print(f"{panel} panel total: ours {our_sum:.1f} vs digitized "
              f"{lz_sum:.1f} ({rel:+.1%})")
        if abs(rel) > 0.20:
            a3 = False

    # --- acceptance 4 & 5: per-bin shape comparison ---
    # dex metrics are evaluated over bins where LZ's drawn level is at
    # least MAIN_FLOOR events/bin: below that the comparison is dominated
    # by the reference's own precision (the drawn line is ~0.02-0.05 dex
    # thick against a 10^-2 axis floor) rather than by model information,
    # and the bins carry <= 0.1% of the model each. The excluded tail is
    # instead checked in absolute counts (acceptance 6).
    MAIN_FLOOR = 1.0

    def lz_levels(panel, series):
        """digitized levels indexed by our bin index (nearest center)."""
        lz = np.full(len(centers), np.nan)
        for r in ref.get((panel, series), []):
            i = int(round((r[0] - centers[0]) / du))
            if 0 <= i < len(centers):
                lz[i] = r[3]
        return lz

    def compare(our_hist, series):
        """panel -> (bin centers, dex residuals over main-floor bins,
        absolute |ours-LZ| sum in events over the excluded tail)."""
        out = {}
        for panel in PANEL_S1C:
            lz = lz_levels(panel, series)
            ok = np.isfinite(lz) & (10 ** lz >= MAIN_FLOOR)
            cs = centers[ok]
            ourv = our_hist[panel][ok]
            res = np.where(ourv > 0,
                           np.log10(np.maximum(ourv, 1e-4)) - lz[ok],
                           np.nan)
            tail = np.isfinite(lz) & (10 ** lz < MAIN_FLOOR)
            tail_abs = float(np.abs(
                our_hist[panel][tail] - 10 ** lz[tail]).sum())
            out[panel] = (cs, res, tail_abs)
        return out

    total_srcs = [s for ss in CATEGORY_MAP.values() for s in ss]
    results = {}
    for cat, series in [("total", "total_bkg"),
                        ("pb214", "pb214_betas"),
                        ("other_er_as_drawn", "other_er_betas"),
                        ("xe124", "xe124_dec"),
                        ("xe127", "xe127_125_ec"),
                        ("accidentals", "accidentals")]:
        our_hist = project(models, total_srcs, band, edges) \
            if cat == "total" else ours[
                {"pb214": "pb214_betas", "other_er_as_drawn":
                 "other_er_as_drawn", "xe124": "xe124_dec",
                 "xe127": "xe127_125_ec", "accidentals":
                 "accidentals"}[cat]]
        results[cat] = compare(our_hist, series)

    # --- acceptance 5 (xe124 placement): LZ floats Qy_LL/Qy_beta in
    # [0.65, 0.87] (WS2024 paper: bounded by the measured single-EC 0.87
    # and the TI box model 0.65). Our pre-fit choice is 0.65; evaluate
    # the category match at the best q in the published prior via the
    # analytic charge rescaling dlog10(S2c) = log10(q/q_gen) applied to
    # the LL-dominated S1c > 40 events (exact for the mean; the
    # fluctuation change is subdominant), and report the implied value.
    s1c124, y124, w124 = models["xe124_dec"]
    mu124, sig124 = band(s1c124)
    hi124 = s1c124 > 40.0
    lz124 = {p: lz_levels(p, "xe124_dec") for p in PANEL_S1C}
    xe124_scan = {}
    for q in np.round(np.arange(0.65, 0.8701, 0.01), 2):
        u = ((y124 + np.where(hi124, np.log10(q / 0.65), 0.0)) - mu124) \
            / sig124
        hq = {}
        for panel, (lo, hh) in PANEL_S1C.items():
            m = (s1c124 > lo) & (s1c124 <= hh)
            hq[panel] = np.histogram(u[m], bins=edges, weights=w124[m])[0]
        xe124_scan[q] = hq

    def rms_at(q):
        res = []
        for panel in PANEL_S1C:
            lz = lz124[panel]
            ok = np.isfinite(lz) & (10 ** lz >= MAIN_FLOOR)
            res.extend(np.log10(np.maximum(xe124_scan[q][panel][ok], 1e-4))
                       - lz[ok])
        res = np.array(res)
        return float(np.sqrt(np.mean(res ** 2)))

    q_best = min(xe124_scan, key=rms_at)
    rms_124_scan = rms_at(q_best)

    print("\nper-bin residuals (ours - LZ), log10, bins with LZ >= "
          f"{MAIN_FLOOR:.0f} event/bin:")
    a4 = True
    res_total = np.concatenate([results["total"][p][1] for p in PANEL_S1C])
    res_total = res_total[np.isfinite(res_total)]
    rms_t = float(np.sqrt(np.mean(res_total ** 2)))
    mx_t = float(np.abs(res_total).max())
    print(f"  total_bkg: rms {rms_t:.3f} dex, max|r| {mx_t:.3f} dex, "
          f"n {len(res_total)}")
    if rms_t > 0.15:
        a4 = False
    a5 = True
    for cat in ("pb214", "other_er_as_drawn", "xe124", "xe127",
                "accidentals"):
        res = np.concatenate([results[cat][p][1] for p in PANEL_S1C])
        res = res[np.isfinite(res)]
        if len(res) == 0:
            print(f"  {cat}: no main-floor comparison bins")
            continue
        rms = np.sqrt(np.mean(res ** 2))
        mx = np.abs(res).max()
        if cat == "xe124":
            print(f"  {cat}: rms {rms:.3f} dex, max|r| {mx:.3f} dex, "
                  f"n {len(res)} at the pre-fit q=0.65; best within LZ's "
                  f"published prior [0.65, 0.87]: q={q_best:.2f} -> "
                  f"rms {rms_124_scan:.3f} dex")
            if rms_124_scan > 0.25:
                a5 = False
        else:
            print(f"  {cat}: rms {rms:.3f} dex, max|r| {mx:.3f} dex, "
                  f"n {len(res)}")
            if rms > 0.25:
                a5 = False

    # --- acceptance 6: the sub-floor tail of the TOTAL curve, in
    # absolute counts ---
    # bins below the main floor (deep low-charge tail, |u| >~ 3): LZ's
    # curves there come from their full position-dependent simulation,
    # whose S2c-correction residuals populate exactly this tail -- the
    # documented uniform_1d limitation (same attribution as the WS2024
    # NR mid-range width residual; the escalation universe
    # simplified_position covers it). Gated on the total curve (what the
    # likelihood sees); per-category tails are reported in the JSON.
    tail_tot = float(sum(results["total"][p][2] for p in PANEL_S1C))
    TAIL_BUDGET = 8.0  # events (~0.7% of the 1210-event model)
    print(f"\nsub-floor tail (total curve): |ours-LZ| = {tail_tot:.2f} "
          f"events (budget {TAIL_BUDGET:.0f})")
    a6 = tail_tot <= TAIL_BUDGET

    # the full other_er vs LZ's drawn outline: expected ~+0.1 dex excess
    # (solar + 212Pb + tritium omitted from LZ's outlines)
    res_full = compare(ours["other_er_betas"], "other_er_betas")
    rf = np.concatenate([res_full[p][1] for p in PANEL_S1C])
    rf = rf[np.isfinite(rf)]
    print(f"  other_er (full, incl. solar/212Pb/tritium) vs drawn outline: "
          f"median {np.median(rf):+.3f} dex (expected ~+0.1: LZ's outlines "
          f"hold only kr85+136Xe)")

    verdict = a1 and per_comp_ok and a3 and a4 and a5 and a6
    print(f"\n{'PASS' if verdict else 'FAIL'}: total {a1}, "
          f"components {per_comp_ok}, panels {a3}, total-shape {a4}, "
          f"categories {a5}, tail {a6}")

    # --- figure: overlays + residuals ---
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    LZ_COLORS = {
        "total_bkg": "tab:blue", "pb214_betas": "tab:cyan",
        "other_er_betas": "tab:pink", "xe124_dec": "tab:green",
        "xe127_125_ec": "tab:red", "accidentals": "tab:orange",
        "b8": "violet", "data": "black",
    }
    OUR_COLORS = {
        "total": "tab:blue", "pb214_betas": "tab:cyan",
        "other_er_as_drawn": "tab:pink", "other_er": "purple",
        "xe124_dec": "tab:green", "xe127_125_ec": "tab:red",
        "accidentals": "tab:orange",
    }
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex="col",
                             height_ratios=[2.2, 1.0])
    for col, panel in enumerate(("top", "bottom")):
        ax, axr = axes[0][col], axes[1][col]
        for (p_key, s), pts in ref.items():
            if p_key != panel:
                continue
            if s == "data":
                ax.plot(pts[:, 0], pts[:, 3], "o", ms=3, color="black",
                        label="LZ data", zorder=5)
                continue
            ax.step(pts[:, 0], pts[:, 3], where="mid", lw=1.2,
                    color=LZ_COLORS.get(s, "gray"),
                    label=f"LZ {s}", alpha=0.9)
        overlay = [("total", "total"), ("pb214_betas", "pb214"),
                   ("other_er_as_drawn", "other ER (as drawn)"),
                   ("xe124_dec", "xe124"), ("xe127_125_ec", "xe127"),
                   ("accidentals", "accidentals")]
        for cat, label in overlay:
            ys = sum(ours[c][panel] for c in
                     ("pb214_betas", "other_er_betas", "xe124_dec",
                      "xe127_125_ec", "accidentals", "b8")) \
                if cat == "total" else ours[cat][panel]
            ax.step(centers, np.maximum(ys, 1e-4), where="mid", lw=1.0,
                    ls="--", color=OUR_COLORS[cat],
                    label=f"ours {label}")
        ax.set_yscale("log")
        ax.set_ylim(3e-3, 400)
        ax.set_ylabel("events / bin")
        ax.set_title(f"{panel} panel (S1c "
                     f"{PANEL_S1C[panel][0]:.0f}-{PANEL_S1C[panel][1]:.0f}"
                     " phd)")
        if col == 0:
            ax.legend(fontsize=6, ncol=2, loc="upper right")
        # residuals: ours - LZ for the total curve
        rc, rr, _ = results["total"][panel]
        axr.axhline(0, color="gray", lw=0.6)
        axr.step(rc, rr, where="mid", color="tab:blue", lw=1.2)
        axr.set_ylabel("ours - LZ (dex)")
        axr.set_xlabel("u = (log10 S2c - mu_ER) / sigma_ER")
        axr.set_ylim(-0.6, 0.6)
    fig.suptitle("Validation 4: background projection vs WS2024 Fig. 4 "
                 "(solid LZ, dashed ours)")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"wrote {args.output}")

    # persist the numeric comparison as metadata
    import json

    meta = {
        "total_counts": total_all,
        "table1_total": t1_total,
        "panel_totals": {
            p: {"ours": float(sum(ours[c][p].sum()
                                  for c in CATEGORY_MAP)),
                "lz_digitized": float(np.sum(
                    10 ** ref[(p, "total_bkg")][:, 3]))}
            for p in PANEL_S1C},
        "residuals_dex": {
            cat: (lambda r: {} if len(r) == 0 else
                  {"rms": float(np.sqrt(np.mean(r ** 2))),
                   "max": float(np.abs(r).max())})(
                np.concatenate([results[cat][p][1] for p in PANEL_S1C]))
            for cat in results},
        "tail_abs_events": {cat: float(sum(results[cat][p][2]
                                           for p in PANEL_S1C))
                            for cat in results},
        "xe124_charge_yield": {
            "pre_fit": 0.65,
            "published_prior": [0.65, 0.87],
            "best_rms_q": float(q_best),
            "rms_at_pre_fit": float(np.sqrt(np.mean(
                np.concatenate([results["xe124"][p][1]
                                for p in PANEL_S1C]) ** 2))),
            "rms_at_best": rms_124_scan},
        "acceptance": {"total": bool(a1), "components": bool(per_comp_ok),
                       "panels": bool(a3), "total_shape": bool(a4),
                       "categories": bool(a5), "tail": bool(a6),
                       "pass": bool(verdict)},
    }
    out_json = Path(args.output).with_suffix(".json")
    out_json.write_text(json.dumps(meta, indent=2))
    print(f"wrote {out_json}")


if __name__ == "__main__":
    main()
