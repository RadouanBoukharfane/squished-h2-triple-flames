"""Fits sT/sT,th = A - b eps* per condition and transport, with eps* measured on the refined cold runs (COLDR-*) and
with eps* reconstructed (same points), from eps_refined_COLDR.json (this directory) and ../data/results_current.json.
Only series whose refined cold runs are all available are fitted. Usage: python fit_eps_measured.py"""
import csv, json, os
import numpy as np

W = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(W, "..", "data", "results_current.json")))
E = {r["cold"]: r for r in json.load(open(os.path.join(W, "eps_refined_COLDR.json"))) if "error" not in r}
# sT,th of all series read from results_current.json (MA harmonised on Cantera, 30/10/2026 14h10)
# paper series -> {point: refined cold run}
MAP = {
    "ma_1bar": {"C1": "C1", "C4": "C4", "C5": "C5", "C6": "C6"},
    "ma_10bar": {"D20": "p10-D20", "D5b": "p10-D5b", "D5": "p10-D5", "D3b": "p10-D3b"},
    "ma_20bar": {"D20": "p20-D20", "D7": "p20-D7", "D5": "p20-D5", "D3": "p20-D3"},
    "ma_1bar579K": {"D20": "B3-MA-D20", "D5b": "B3-MA-D5b", "D5": "B3-MA-D5"},
    "sc_1bar": {"C1": "C1-SC", "C3": "C3-SC", "C4": "C4-SC", "A6-tanh": "A6-tanh-SC", "A6-erf": "A6-erf-SC",
                  "C5": "C5-SC", "C6": "C6-SC"},
    "sc_10bar": {k: "B1-p10-%s-SC" % k for k in ("D20", "D7", "D5", "D4", "D3")},
    "sc_20bar": {k: "B1-p20-%s-SC" % k for k in ("D20", "D7", "D5", "D4", "D3")},
    "sc_1bar579K": {"D32": "B3-D32-SC", "D7": "B3-D7-SC", "D3": "B3-D3-SC"},
    "sc_10bar300K": {"D32": "B2-D32-SC", "D7": "B2-D7-SC", "D5": "B2-D5-SC"},
    "ma_10bar300K": {k: "B2-MA-" + k for k in ("D20", "D5b", "D5")},
    "ma_5bar": {k: "B5-MA-" + k for k in ("D20", "D5b", "D5")},
}


def fit(e, s):
    (bb, aa), cov = np.polyfit(e, s, 1, cov=True)
    return aa, -bb, np.sqrt(cov[0, 0]), np.sqrt(np.mean((s - aa - bb * e) ** 2))


res = {}
for ser, m in MAP.items():
    th = D["sTth"][ser]
    pts = [(n, "COLDR-" + c) for n, c in m.items() if D[ser].get(n, {}).get("dplus")]
    miss = [c for n, c in pts if c not in E]
    if miss:
        print("%-15s incomplete, missing: %s" % (ser, " ".join(miss))); continue
    s = np.array([D[ser][n]["sT"] / th for n, c in pts])
    em = np.array([E[c]["eps_measured"] for n, c in pts]); er = np.array([E[c]["eps_reconstructed"] for n, c in pts])
    res[ser] = dict(meas=fit(em, s), rec=fit(er, s), n=len(pts))
    print("%-15s n=%d  measured: A=%.3f b=%.3f±%.3f rms=%.3f | reconstructed: A=%.3f b=%.3f±%.3f rms=%.3f"
          % ((ser, len(pts)) + res[ser]["meas"] + res[ser]["rec"]))
print()
for code in ("ma", "sc"):
    b1 = res.get(code + "_1bar")
    if not b1:
        continue
    for ser in res:
        if ser.startswith(code + "_") and ser != code + "_1bar":
            p = ser.split("_", 1)[1]
            pr = D["prediction_b_ratio"].get(p, {})
            print("%-15s b/b1 measured %.2f   reconstructed %.2f   predicted lambda Ma_c %s, Ma_u %s"
                  % (ser, res[ser]["meas"][1] / b1["meas"][1], res[ser]["rec"][1] / b1["rec"][1],
                     pr.get("Ma_c"), pr.get("Ma_u")))
json.dump({k: dict(n=v["n"], **{t: dict(zip(("A", "b", "sb", "rms"), map(float, v[t]))) for t in ("meas", "rec")})
           for k, v in res.items()}, open(os.path.join(W, "fits_eps_measured.json"), "w"), indent=1)
