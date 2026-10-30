"""eps* measured on the computed cold layer (review of 29/10/2026, point 5), against the reconstructed one.
For each cold run COLDF-<case> (or COLD-<case>): mixture fraction z(x) along x at the tip height of the reactive case
(active_control.height), finest level (fextract), interpolated monotonically (PCHIP) on a fine grid; s(x) = sL(z(x)) with the Cantera polynomial sL(z) of the condition
(same tables as make_paper_figures.py); eps* = sqrt(|s''|/s) at the maximum of s, x in deltaL units, s'' from a
quadratic fit on a window of half-width 0.4 R around the maximum (iterated). The same estimator is applied to the
reconstructed profile (tanh of thickness Delta0+ diffused by a Gaussian to the measured Delta+) for comparison.
Usage (PELE_HOME set, fextract built; run directories under $RUNS):
  python eps_measured.py OUT.json COLDF-C4 COLDF-S6-tanh ...     (run directories under $RUNS)"""
import csv, glob, json, os, re, subprocess, sys
import numpy as np
from scipy.signal import fftconvolve
from scipy.interpolate import PchipInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PELE = os.environ.get("RUNS", os.path.join(ROOT, "cases"))
PH = os.environ.get("PELE_HOME", "")                # only needed to extract profiles (fextract)
FX = os.path.join(PH, "Submodules/PelePhysics/Submodules/amrex/Tools/Plotfile/fextract.gnu.ex")
zst = 0.02830074092068505
phi2z = lambda phi: phi * zst / (1.0 - zst + phi * zst)
TABLE = {(1, 300): ("p1_T300", 25.0), (10, 579): ("p10_T579", 2.0), (20, 706): ("p20_T706", 0.92),
         (1, 579): ("p1_T579", 34.0), (10, 300): ("p10_T300", 1.6),
         (5, 475): ("p5_T475", 4.1)}


def sL_poly(cond, dx):
    rows = [r for r in csv.DictReader(open(os.path.join(ROOT, "theory", "sL_cantera_phi.csv"))) if r["condition"] == cond]
    d = {float(r["phi"]): float(r["sL_cantera_m_s"]) for r in rows}
    phi = np.array(sorted(d))
    return np.polyfit(phi2z(phi), [d[p] for p in phi], 4)


def eps_fit(x, s):
    """sqrt(|s''|/s) at the maximum of s (x in deltaL units), iterated local quadratic fit."""
    m = np.argmax(s)
    w = 1.0
    for _ in range(8):
        sel = np.abs(x - x[m]) < w
        if sel.sum() < 5:
            sel = np.abs(x - x[m]) <= np.sort(np.abs(x - x[m]))[4]
        c = np.polyfit(x[sel] - x[m], s[sel], 2)
        eps = np.sqrt(abs(2 * c[0]) / c[2])
        w_new = 0.4 / eps
        if abs(w_new - w) < 0.02 * w:
            break
        w = w_new
    return eps, int(sel.sum())


def layer(x, d0p, dp):
    h = x[1] - x[0]
    z = 0.5 * (1 + np.tanh(2 * x / d0p))
    dd = np.sqrt(max(dp ** 2 - d0p ** 2, 0.0))
    if dd > 0:
        sd = dd / np.sqrt(2 * np.pi)
        k = np.exp(-0.5 * (np.arange(-int(8 * sd / h), int(8 * sd / h) + 1) * h / sd) ** 2)
        z = fftconvolve(z, k / k.sum(), mode="same")
    return z


def one(cold):
    case = re.sub(r"^COLD[FR]?-", "", cold)
    if not os.path.isdir(os.path.join(PELE, case)) and os.path.isdir(os.path.join(PELE, "B1-" + case)):
        case = "B1-" + case  # cold runs COLDF-p10-* / COLDF-p20-* of the reactive cases B1-p10-* / B1-p20-*
    inp = open(os.path.join(PELE, case, "input.case" if os.path.exists(os.path.join(PELE, case, "input.case"))
                            else "input.c1")).read()
    g = lambda k: re.search(r"^%s\s*=\s*(\S+)" % re.escape(k), inp, re.M).group(1)
    W = float(re.search(r"^geometry.prob_hi\s*=\s*(\S+)", inp, re.M).group(1))
    dL = 317e-6 * W / 0.075
    P = float(g("prob.P_mean")); p = 1 if P < 2e5 else int(round(P / 1e5))
    tu = int(round(float(g("prob.T_in"))))
    poly = sL_poly(*TABLE[(p, tu)])
    h = float(g("active_control.height"))
    d0p = float(g("prob.delta0")) / dL
    pf = sorted(glob.glob(os.path.join(PELE, cold, "plt?????")))[-1]
    out = os.path.join(PELE, cold, "slice_tip.dat")
    subprocess.run([FX, "-d", "0", "-y", "%.8e" % h, "-v", "mixture_fraction", "-s", out, pf], check=True,
                   stdout=subprocess.DEVNULL)
    x, z = np.loadtxt(out, unpack=True)
    xd = x / dL
    dp = 1.0 / np.max(np.gradient(z, xd))
    # z(x) on the cold-run grid (about 0.25 deltaL at the tip height) interpolated monotonically (PCHIP) on a fine grid
    xd, iu = np.unique(xd, return_index=True); z = z[iu]
    xf = np.arange(xd[0], xd[-1], 0.002)
    zf = PchipInterpolator(xd, z)(xf)
    ok = (zf > 0.005) & (zf < phi2z(2.0))
    eps_m, n_m = eps_fit(xf[ok], np.polyval(poly, zf[ok]))
    xr = np.linspace(-80, 80, 160001)
    zr = layer(xr, d0p, dp)
    okr = (zr > 0.005) & (zr < phi2z(2.0))
    eps_r, _ = eps_fit(xr[okr], np.polyval(poly, zr[okr]))
    dx_fine = np.min(np.diff(xd))
    return dict(cold=cold, case=case, p=p, Tu=tu, d0p=d0p, dplus=dp, eps_measured=eps_m, eps_reconstructed=eps_r,
                ratio=eps_m / eps_r, points_in_fit=n_m, dx_over_dL=dx_fine, Yplus=h / dL)


if __name__ == "__main__":
    res = []
    for c in sys.argv[2:]:
        try:
            r = one(c)
            print("%-24s p=%2d Tu=%3d Y+=%5.1f d0p=%5.2f D+=%5.2f eps_meas=%.3f eps_rec=%.3f ratio=%.3f (n=%d, dx=%.3f dL)"
                  % (r["cold"], r["p"], r["Tu"], r["Yplus"], r["d0p"], r["dplus"], r["eps_measured"],
                     r["eps_reconstructed"], r["ratio"], r["points_in_fit"], r["dx_over_dL"]), flush=True)
        except Exception as e:                          # noqa: BLE001
            r = dict(cold=c, error=str(e)); print(c, "ERROR", e, flush=True)
        res.append(r)
    json.dump(res, open(sys.argv[1], "w"), indent=1)
