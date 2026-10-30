"""Post-processing of tip_analysis.json (../tools/tip_analysis.py), review of 29/10/2026, points 1 and 5.
(1) tip Markstein number per series: deficit 1 - s_c/sL against the Karlovitz number of the total stretch at the
    unburnt edge, Ka = K deltaL/sL (sL = s_m, maximum of sL(z)), fit through the origin and with intercept; compared
    with Ma_c,delta of the counterflow flames (markstein_defs.json, MA; constant Sc: prediction_b/lambda). Two-Markstein
    fit (strain and curvature parts separately) for information.
(2) sT/sTth = A - b eps*, with eps* measured in the reacting flow k deltaL below the tip (k = 2, 4, 8): b, rms,
    b/b1, b/(lambda Ma_c), ratio MA/SC; compared with the reconstructed eps* (fit_fixed_A.json, free A).
Usage: python tip_post.py [tip_analysis.json]"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
from eps_measured import TABLE, sL_poly                                    # noqa: E402

W = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(W))
T = {r["case"]: r for r in json.load(open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(W, "tip_analysis.json")))
     if "error" not in r}
FA = json.load(open(os.path.join(W, "fits_fixed_A.json")))
D = json.load(open(os.path.join(W, "..", "data", "results_current.json")))
MD = {(r["p"], r["Tu"]): r for r in json.load(open(os.path.join(ROOT, "theory", "markstein_defs.json")))
      if r.get("mech") == "sandiego"}
MD[(5, 475)] = [r for r in json.load(open(os.path.join(ROOT, "theory", "markstein_sweep.json"))) if r["p"] == 5][0]
COND = {"1bar": (1, 300), "10bar": (10, 579), "20bar": (20, 706), "1bar579K": (1, 579), "10bar300K": (10, 300),
        "5bar": (5, 475)}
CASES = {"ma_1bar": ["C1", "C2", "C4", "C5", "C6"],
         "ma_10bar": ["B1-p10-" + k for k in ("D20", "D5b", "D5", "D3b")],
         "ma_20bar": ["B1-p20-" + k for k in ("D20", "D7", "D5", "D3")],
         "ma_1bar579K": ["B3-MA-" + k for k in ("D20", "D5b", "D5")],
         "ma_10bar300K": ["B2-MA-" + k for k in ("D20", "D5b", "D5")],
         "ma_5bar": ["B5-MA-" + k for k in ("D20", "D5b", "D5")],
         "sc_1bar": ["C1-SC", "C3-SC", "C4-SC", "A6-tanh-SC", "A6-erf-SC", "C5-SC", "C6-SC"],
         "sc_10bar": ["B1-p10-%s-SC" % k for k in ("D20", "D7", "D5", "D4", "D3")],
         "sc_20bar": ["B1-p20-%s-SC" % k for k in ("D20", "D7", "D5", "D4", "D3")],
         "sc_1bar579K": ["B3-%s-SC" % k for k in ("D32", "D7", "D3")],
         "sc_10bar300K": ["B2-%s-SC" % k for k in ("D32", "D7", "D5")]}
th = D["sTth"]                       # MA harmonised on Cantera in results_current.json (30/10/2026 14h10)
# reference speed: s_m = max of the Cantera sL(z) of the condition (the normalisation of the paper); the vertical
# through the tip does not give a usable local z (preferential diffusion near the flame, layer deflected upstream)
SM = {}
for c, r in T.items():
    zz = np.linspace(0.005, 0.07, 2000)
    SM[c] = float(np.polyval(sL_poly(*TABLE[(r["p"], int(round(r["Tu"])))]), zz).max())
    r["sc_over_sm"] = r["s_c"] / SM[c]
    for q in ("Ka_u", "Ka_s_u", "Ka_k_u", "Ka_1000"):
        r[q + "_m"] = r[q] * r["sL"] / SM[c]
V = {}                                                     # final inlet velocity = sT
for c in [c for v in CASES.values() for c in v]:
    rows = [l.split() for l in open(os.path.join(os.environ.get("RUNS", os.path.join(ROOT, "cases")), c, "AC_History.dat")) if l.strip()]
    V[c] = float(rows[-1][2])


def ma_c(ser):
    code, p = ser.split("_", 1)
    m = MD[COND[p]]
    if code == "ma":
        return m["Ma_c"], m["lam"]
    pb = D["prediction_b"].get(ser)
    return (pb / m["lam"] if pb else np.nan), m["lam"]


print("(1) tip Markstein number: deficit 1 - s_c/sL against Ka of the total stretch at the unburnt edge")
print("%-16s %-3s %-9s %-14s %-9s %-26s %s" % ("series", "n", "Ma_tip0", "Ma_tip (A)", "Ma_c", "Ma_tip0/Ma_c", "two-Markstein M_s, M_k"))
out = {}
for ser, cs in CASES.items():
    cs = [c for c in cs if c in T]
    y = np.array([1 - T[c]["sc_over_sm"] for c in cs])
    x = np.array([T[c]["Ka_u_m"] for c in cs])
    xs = np.array([T[c]["Ka_s_u_m"] for c in cs]); xk = np.array([T[c]["Ka_k_u_m"] for c in cs])
    m0 = float(np.sum(x * y) / np.sum(x * x))
    m1, a1 = np.polyfit(x, y, 1)
    M2 = np.linalg.lstsq(np.c_[xs, xk], y, rcond=None)[0]
    mc, lam = ma_c(ser)
    out[ser] = dict(cases=cs, Ma_tip0=m0, Ma_tip=float(m1), intercept=float(a1), Ma_c=mc, lam=lam,
                    M_s=float(M2[0]), M_k=float(M2[1]),
                    points=[dict(case=c, Ka=float(a), deficit=float(b), Ka_s=float(s), Ka_k=float(k), phi_u=T[c]["phi_u"])
                            for c, a, b, s, k in zip(cs, x, y, xs, xk)])
    print("%-16s %-3d %-9.3f %.3f (%+.3f) %-9.3f %-26.2f %.3f, %.3f" % (ser, len(cs), m0, m1, a1, mc, m0 / mc, M2[0], M2[1]))
    for c, a, b in zip(cs, x, y):
        print("      %-15s Ka=%.3f (s %.3f, k %.3f)  1-sc/sL=%+.3f  phi_u=%.2f  Ka_1000=%.3f  kappa dL=%.3f" %
              (c, a, T[c]["Ka_s_u_m"], T[c]["Ka_k_u_m"], b, T[c]["phi_u"], T[c]["Ka_1000_m"], T[c]["kappa_dL_1000"]))

print("\n(2) sT/sTth = A - b eps*, eps* measured in the reacting flow k deltaL below the tip")
fits = {}
for k in ("rec", 2, 4, 8):
    for ser, cs in CASES.items():
        if k == "rec":
            f = FA[ser]["free"]; fits[(ser, k)] = (f["A"], f["b"], f["rms"]); continue
        cs = [c for c in cs if c in T]
        e = np.array([T[c]["eps_up%d" % k] for c in cs]); s = np.array([V[c] / th[ser] for c in cs])
        bb, aa = np.polyfit(e, s, 1)
        fits[(ser, k)] = (float(aa), float(-bb), float(np.sqrt(np.mean((s - aa - bb * e) ** 2))))
for k in ("rec", 2, 4, 8):
    print("-- eps* %s" % ("reconstructed (paper)" if k == "rec" else "reacting flow, %d deltaL below the tip" % k))
    for ser in CASES:
        A, b, rms = fits[(ser, k)]
        b1 = fits[(ser.split("_")[0] + "_1bar", k)][1]
        mc, lam = ma_c(ser)
        other = ("sc_" if ser.startswith("ma") else "ma_") + ser.split("_", 1)[1]
        rt = "MA/SC %.2f" % (b / fits[(other, k)][1]) if ser.startswith("ma") and (other, k) in fits else ""
        print("   %-16s A=%.3f b=%.3f rms=%.4f  b/b1=%.2f  b/(lam Ma_c)=%.2f  %s" % (ser, A, b, rms, b / b1, b / (lam * mc), rt))
print("predicted b/b1 (Ma_c, Ma_u):", D["prediction_b_ratio"], "; lambda Ma_c MA/SC = %.2f (1 atm)"
      % (ma_c("ma_1bar")[0] / ma_c("sc_1bar")[0]))
json.dump(dict(markstein=out, eps_fits={"%s|%s" % key: v for key, v in fits.items()}),
          open(os.path.join(W, "tip_post.json"), "w"), indent=1)
