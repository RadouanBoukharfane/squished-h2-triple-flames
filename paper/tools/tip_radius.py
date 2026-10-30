"""Tip radius of the PeleLMeX triple flames (figure 4): radius R of the T = 1000 K isotherm at the leading edge
(lowest point), from a parabola y = y0 + c (x - x0)^2 + ... fitted on the isotherm near its vertex, on the finest
AMR level. The fit window is iterated to |x - x0| < 0.5 R. Also stores the isotherm around the tip in units of
deltaL for the shape panel. deltaL = 317 um x ls, ls = (domain width)/(75 mm).
Usage (compute node): python tip_radius.py OUT.json CASE [CASE ...]   (CASE = run directory under $RUNS)"""
import glob, json, os, sys
import numpy as np
import contourpy
import yt

yt.set_log_level(50)
PELE = os.environ.get("RUNS", os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "cases"))
T_TIP = float(os.environ.get("T_TIP", "1000"))   # isotherm (env T_TIP)


def grid(ds, lev, lo, hi):
    dx = float(ds.domain_width[0]) / (ds.domain_dimensions[0] * 2 ** lev)
    i0, j0 = max(int(lo[0] / dx), 0), max(int(lo[1] / dx), 0)
    ni = min(int((hi[0] - lo[0]) / dx), ds.domain_dimensions[0] * 2 ** lev - i0)
    nj = min(int((hi[1] - lo[1]) / dx), ds.domain_dimensions[1] * 2 ** lev - j0)
    cg = ds.covering_grid(lev, left_edge=[i0 * dx, j0 * dx, 0.0], dims=[ni, nj, 1])
    T = np.array(cg[("boxlib", "temp")])[:, :, 0].T
    x = (i0 + 0.5 + np.arange(ni)) * dx
    y = (j0 + 0.5 + np.arange(nj)) * dx
    return x, y, T, dx


def tip(case):
    pf = sorted(glob.glob(os.path.join(PELE, case, "plt?????")))[-1]
    ds = yt.load(pf)
    W, H = float(ds.domain_width[0]), float(ds.domain_width[1])
    ls = W / 0.075
    dL = 317e-6 * ls
    # coarse search of the leading edge (level 1, whole domain)
    x, y, T, _ = grid(ds, 1, (0.0, 0.0), (W, H))
    hot = T > T_TIP
    j, i = np.unravel_index(np.argmin(np.where(hot, y[:, None] + 0 * x[None, :], 1e9)), T.shape)
    xt, yt_ = x[i], y[j]
    # finest level around the tip
    lev = ds.max_level
    x, y, T, dx = grid(ds, lev, (xt - 8 * dL, yt_ - 3 * dL), (xt + 8 * dL, yt_ + 10 * dL))
    lines = contourpy.contour_generator(x, y, T).lines(T_TIP)
    pts = np.concatenate([l for l in lines if len(l) > 10])
    k = np.argmin(pts[:, 1])
    x0, y0 = pts[k]
    w = 1.0 * dL
    for _ in range(6):
        sel = np.abs(pts[:, 0] - x0) < w
        c = np.polyfit(pts[sel, 0] - x0, pts[sel, 1], 4)
        xv = np.linspace(-w, w, 2001)
        yv = np.polyval(c, xv)
        kv = np.argmin(yv)
        x0, y0 = x0 + xv[kv], yv[kv]
        curv = np.polyval(np.polyder(c, 2), xv[kv])           # y'' at the vertex (y' = 0): kappa = y''
        R = 1.0 / curv
        w_new = max(0.5 * R, 4 * dx)
        if abs(w_new - w) < 0.02 * w:
            break
        w = w_new
    sel = (np.abs(pts[:, 0] - x0) < 6 * dL) & (pts[:, 1] - y0 < 8 * dL)
    shp = pts[sel]
    shp = shp[np.argsort(shp[:, 0])][:: max(1, sel.sum() // 400)]
    return dict(case=case, plt=os.path.basename(pf), ls=ls, dL=dL, dx=dx, R=float(R), R_over_dL=float(R / dL),
                eps_tip=float(dL / R), x_tip=float(x0), y_tip=float(y0), window=float(w),
                shape=[[float((a - x0) / dL), float((b - y0) / dL)] for a, b in shp])


if __name__ == "__main__":
    out = []
    for case in sys.argv[2:]:
        try:
            r = tip(case)
            print("%-14s %s  dL = %.2f um  dx = %.2f um  R = %.1f um  R/dL = %.3f  dL/R = %.3f  window %.2f dL"
                  % (case, r["plt"], r["dL"] * 1e6, r["dx"] * 1e6, r["R"] * 1e6, r["R_over_dL"], r["eps_tip"],
                     r["window"] / r["dL"]), flush=True)
        except Exception as e:                       # noqa: BLE001
            r = dict(case=case, error=str(e))
            print(case, "ERROR", e, flush=True)
        out.append(r)
    json.dump(out, open(sys.argv[1], "w"), indent=1)
