"""Check of tools/tip_radius.py (28/10/2026): for each case and isotherm, tip radius from (i) the iterated polynomial
fit with windows |x - x0| < f R, f = 0.25, 0.5, 1, (ii) a least-squares circle (Kasa) on the isotherm points within
an arc of +-0.5 R around the vertex; resolution R/dx; vertex position relative to the stoichiometric line and to
the layer centre; number of separate isotherm pieces below the vertex height + 2 deltaL.
Usage (compute node): python tip_radius_check.py OUT.json T1,T2,... CASE [CASE ...]"""
import glob, json, os, sys
import numpy as np
import contourpy
import yt

yt.set_log_level(50)
PELE = os.environ.get("RUNS", os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "cases"))
ZST = 0.02830074092068505


def grid(ds, lev, lo, hi, fields):
    dx = float(ds.domain_width[0]) / (ds.domain_dimensions[0] * 2 ** lev)
    i0, j0 = max(int(lo[0] / dx), 0), max(int(lo[1] / dx), 0)
    ni = min(int((hi[0] - lo[0]) / dx), ds.domain_dimensions[0] * 2 ** lev - i0)
    nj = min(int((hi[1] - lo[1]) / dx), ds.domain_dimensions[1] * 2 ** lev - j0)
    cg = ds.covering_grid(lev, left_edge=[i0 * dx, j0 * dx, 0.0], dims=[ni, nj, 1])
    out = [np.array(cg[("boxlib", f)])[:, :, 0].T for f in fields]
    return (i0 + 0.5 + np.arange(ni)) * dx, (j0 + 0.5 + np.arange(nj)) * dx, out, dx


def poly_R(pts, x0, w0, frac):
    w = w0
    for _ in range(10):
        sel = np.abs(pts[:, 0] - x0) < w
        if sel.sum() < 6:
            return np.nan, x0
        c = np.polyfit(pts[sel, 0] - x0, pts[sel, 1], 4)
        xv = np.linspace(-w, w, 4001)
        kv = np.argmin(np.polyval(c, xv))
        x0 = x0 + xv[kv]
        R = 1.0 / np.polyval(np.polyder(c, 2), xv[kv])
        w_new = frac * R
        if abs(w_new - w) < 0.01 * w:
            break
        w = w_new
    return R, x0


def circle_R(pts, x0, y0, R):
    d = np.hypot(pts[:, 0] - x0, pts[:, 1] - y0)
    sel = d < 0.5 * R
    x, y = pts[sel, 0], pts[sel, 1]
    A = np.c_[2 * x, 2 * y, np.ones_like(x)]
    sol = np.linalg.lstsq(A, x ** 2 + y ** 2, rcond=None)[0]
    return np.sqrt(sol[2] + sol[0] ** 2 + sol[1] ** 2), sel.sum()


def check(case, temps):
    pf = sorted(glob.glob(os.path.join(PELE, case, "plt?????")))[-1]
    ds = yt.load(pf)
    W, H = float(ds.domain_width[0]), float(ds.domain_width[1])
    dL = 317e-6 * W / 0.075
    x, y, (T,), _ = grid(ds, 1, (0.0, 0.0), (W, H), ["temp"])
    j, i = np.unravel_index(np.argmin(np.where(T > 1000, y[:, None] + 0 * x[None, :], 1e9)), T.shape)
    xt, yt_ = x[i], y[j]
    lev = ds.max_level
    x, y, (T, Z), dx = grid(ds, lev, (xt - 8 * dL, yt_ - 3 * dL), (xt + 8 * dL, yt_ + 10 * dL), ["temp", "mixture_fraction"])
    res = dict(case=case, dL=dL, dx=dx)
    for Tiso in temps:
        lines = [l for l in contourpy.contour_generator(x, y, T).lines(Tiso) if len(l) > 10]
        pts = np.concatenate(lines)
        k = np.argmin(pts[:, 1]); x0, y0 = pts[k]
        npieces = sum(1 for l in lines if l[:, 1].min() < y0 + 2 * dL)
        Rs = {}
        for f in (0.25, 0.5, 1.0):
            Rs[f], xv = poly_R(pts, x0, 1.0 * dL, f)
        Rc, nc = circle_R(pts, xv, y0 + Rs[0.5], Rs[0.5])
        # mixture fraction at the vertex (bilinear on the grid)
        zi = Z[np.argmin(abs(y - y0)), np.argmin(abs(x - xv))]
        # layer centre: maximum gradient of z along the row of the vertex
        row = Z[np.argmin(abs(y - y0))]
        xc = x[np.argmax(np.abs(np.gradient(row, x)))]
        zst_row = x[np.argmin(abs(row - ZST))]
        res["T%d" % Tiso] = dict(R_025=Rs[0.25] / dL, R_05=Rs[0.5] / dL, R_1=Rs[1.0] / dL, R_circle=Rc / dL,
                                 n_circle=int(nc), R_over_dx=Rs[0.5] / dx, z_vertex=float(zi),
                                 x_vertex_minus_centre=(xv - xc) / dL, x_vertex_minus_xst=(xv - zst_row) / dL,
                                 pieces_near_tip=int(npieces))
    return res


if __name__ == "__main__":
    temps = [float(s) for s in sys.argv[2].split(",")]
    out = []
    for case in sys.argv[3:]:
        try:
            r = check(case, temps)
        except Exception as e:                      # noqa: BLE001
            r = dict(case=case, error=str(e))
        out.append(r)
        print(json.dumps(r), flush=True)
    json.dump(out, open(sys.argv[1], "w"), indent=1)
