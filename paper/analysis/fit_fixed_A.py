"""Review of 29/10/2026, point 4: status of the intercept A of sT/sTth = A - b eps*, and its effect on the slopes.
Same points and same reconstructed eps* as make_paper_figures.py (results_current.json), updated with the runs finished
since (all series now in paper/data/results_current.json).
Fits compared, per transport and condition:
  free   : A and b free (paper);
  A=1    : plateau of (sTth) imposed;
  A=wide : A set to the reduced speed of the widest layer of the series (measured plateau of that series);
  A=com  : one A per transport, common to all conditions (joint least squares, b per condition).
Prints b, b/b1, b/(lambda Ma_c), rms. Writes fits_fixed_A.json. Usage: python fit_fixed_A.py"""
import csv, json, os
import numpy as np
from scipy.signal import fftconvolve

W = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(W))
D = json.load(open(os.path.join(W, "..", "data", "results_current.json")))
zst = 0.02830074092068505
phi2z = lambda phi: phi * zst / (1.0 - zst + phi * zst)


def sL_poly(cond, dx):
    rows = [r for r in csv.DictReader(open(os.path.join(ROOT, "theory", "sL_cantera_phi.csv"))) if r["condition"] == cond]
    d = {float(r["phi"]): float(r["sL_cantera_m_s"]) for r in rows}
    phi = np.array(sorted(d))
    return np.polyfit(phi2z(phi), [d[p] for p in phi], 4)


def layer(x, d0p, dp):
    h = x[1] - x[0]
    z = 0.5 * (1 + np.tanh(2 * x / d0p))
    dd = np.sqrt(max(dp ** 2 - d0p ** 2, 0.0))
    if dd > 0:
        sd = dd / np.sqrt(2 * np.pi)
        k = np.exp(-0.5 * (np.arange(-int(8 * sd / h), int(8 * sd / h) + 1) * h / sd) ** 2)
        z = fftconvolve(z, k / k.sum(), mode="same")
    return z


def eps_star(P, d0p, dp):                     # identical to make_paper_figures.py
    x = np.linspace(-80, 80, 160001)
    h = x[1] - x[0]
    z = layer(x, d0p, dp)
    s = np.polyval(P, z)
    ok = (z > 0.005) & (z < phi2z(2.0))
    m = np.argmax(np.where(ok, s, -1))
    return np.sqrt(abs(np.gradient(np.gradient(s, h), h)[m]) / s[m])


POLY = {"1bar": sL_poly("p1_T300", 25.0), "10bar": sL_poly("p10_T579", 2.0), "20bar": sL_poly("p20_T706", 0.92),
        "1bar579K": sL_poly("p1_T579", 34.0), "10bar300K": sL_poly("p10_T300", 1.6), "5bar": sL_poly("p5_T475", 4.1)}
th = D["sTth"]                       # MA harmonised on Cantera in results_current.json (30/10/2026 14h10)
S = {k: dict(v) for k, v in D.items() if k.startswith(("ma_", "sc_")) and k != "ma_1bar_soret"
     and isinstance(v, dict) and all(isinstance(x, dict) for x in v.values())}
PB = D["prediction_b"]

pts = {}
for ser, m in S.items():
    p = ser.split("_", 1)[1]
    rows = [(n, eps_star(POLY[p], v["d0p"], v["dplus"]), v["sT"] / th[ser], bool(v.get("prov")), v["dplus"])
            for n, v in m.items() if v.get("dplus")]
    pts[ser] = rows


def lsq(e, s, A=None):
    if A is None:
        (bb, aa), _ = np.polyfit(e, s, 1, cov=True) if len(e) > 2 else (np.polyfit(e, s, 1), None)
        return aa, -bb
    b = np.sum((A - s) * e) / np.sum(e * e)
    return A, b


res = {}
for ser, rows in pts.items():
    e = np.array([r[1] for r in rows]); s = np.array([r[2] for r in rows])
    wide = rows[int(np.argmax([r[4] for r in rows]))]
    out = {}
    for tag, A in (("free", None), ("A=1", 1.0), ("A=wide", wide[2])):
        a, b = lsq(e, s, A)
        out[tag] = dict(A=float(a), b=float(b), rms=float(np.sqrt(np.mean((s - a + b * e) ** 2))))
    res[ser] = dict(points=[dict(n=r[0], eps=float(r[1]), red=float(r[2]), prov=r[3]) for r in rows],
                    wide=dict(n=wide[0], red=float(wide[2]), eps=float(wide[1])), **out)
for code in ("ma", "sc"):                  # common A per transport: unknowns A, b_ser
    sers = [k for k in res if k.startswith(code)]
    nrow = sum(len(pts[k]) for k in sers)
    M = np.zeros((nrow, 1 + len(sers))); y = np.zeros(nrow); i = 0
    for j, k in enumerate(sers):
        for r in pts[k]:
            M[i, 0] = 1; M[i, 1 + j] = -r[1]; y[i] = r[2]; i += 1
    sol = np.linalg.lstsq(M, y, rcond=None)[0]
    for j, k in enumerate(sers):
        e = np.array([r[1] for r in pts[k]]); s = np.array([r[2] for r in pts[k]])
        res[k]["A=com"] = dict(A=float(sol[0]), b=float(sol[1 + j]),
                               rms=float(np.sqrt(np.mean((s - sol[0] + sol[1 + j] * e) ** 2))))
TAGS = ("free", "A=1", "A=wide", "A=com")
for code in ("ma", "sc"):
    print("== %s" % {"ma": "mixture-averaged", "sc": "constant Schmidt"}[code])
    print("%-16s %-7s" % ("series", "wide") + "".join("| %-28s" % t for t in TAGS))
    b1 = res[code + "_1bar"]
    for ser in [k for k in res if k.startswith(code)]:
        r = res[ser]; c = []
        for t in TAGS:
            b = r[t]["b"]; pr = PB.get(ser)
            c.append("A=%.3f b=%.3f r=%.2f p=%s rms=%.3f" % (r[t]["A"], b, b / b1[t]["b"],
                                                              "%.2f" % (b / pr) if pr else "  - ", r[t]["rms"]))
        prov = "*" if any(p["prov"] for p in r["points"]) else " "
        print("%-15s%s %.3f  " % (ser, prov, r["wide"]["red"]) + " ".join("| " + x for x in c))
print("r = b/b1 ; p = b/(lambda Ma_c) ; * = provisional points ; predicted b/b1 (Ma_c, Ma_u):", D["prediction_b_ratio"])
json.dump(res, open(os.path.join(W, "fits_fixed_A.json"), "w"), indent=1)
