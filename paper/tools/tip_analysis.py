"""Leading-edge analyses of the steady PeleLMeX triple flames (review of 29/10/2026, points 1 and 5), last plotfile.
Tip = lowest point of the T = 1000 K isotherm (as tip_radius.py). On the finest level around the tip:
 (1) tip Markstein number, same definitions as the counterflow Ma_c of table 1 (theory/markstein_defs.py):
     consumption speed s_c = -int w_O2 dy / (rho_u Y_O2,u) along the vertical through the tip (the normal there), from
     the unburnt side to the first minimum of -w_O2 past its peak (the premixed branch only, not the trailing
     diffusion flame); sL = Cantera polynomial sL(z) at the mixture fraction z_u of the unburnt edge (T = Tu + 10 K)
     on that vertical; stretch K = a_t + s_d kappa with n = -grad T/|grad T| (towards the fresh gas), kappa = div n,
     a_t = div u - n.grad(u).n, s_d = -u.n, evaluated at the isotherms Tu + 10 K (unburnt edge, as the upstream
     strain of the counterflow) and 1000 K; Karlovitz numbers Ka = K deltaL / sL, split into strain and curvature.
 (2) eps* seen by the tip in the reacting flow: z(x) along horizontal lines k deltaL below the tip (k = 2, 4, 8),
     same estimator as eps_measured.py (PCHIP on a fine grid, sL(z), sqrt(|s''|/s) at the maximum of s).
deltaL = 317 um x (domain width)/(75 mm), as the other scripts.
Usage (compute node): python tip_analysis.py OUT.json NPROC CASE [CASE ...]   (CASE = run directory under $RUNS)"""
import glob, json, os, re, subprocess, sys
from multiprocessing import Pool
import numpy as np
import contourpy
from scipy.interpolate import PchipInterpolator, RegularGridInterpolator
import yt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eps_measured import TABLE, FX, sL_poly, eps_fit, phi2z          # noqa: E402

yt.set_log_level(50)
PELE = os.environ.get("RUNS", os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "cases"))
F = ["temp", "x_velocity", "y_velocity", "density", "I_R(O2)", "Y(O2)", "mixture_fraction"]


def cgrid(ds, lev, lo, hi, fields):
    dx = float(ds.domain_width[0]) / (ds.domain_dimensions[0] * 2 ** lev)
    nx, ny = ds.domain_dimensions[0] * 2 ** lev, ds.domain_dimensions[1] * 2 ** lev
    i0, j0 = max(int(lo[0] / dx), 0), max(int(lo[1] / dx), 0)
    ni, nj = min(int((hi[0] - lo[0]) / dx), nx - i0), max(min(int((hi[1] - lo[1]) / dx), ny - j0), 1)
    cg = ds.covering_grid(lev, left_edge=[i0 * dx, j0 * dx, 0.0], dims=[ni, nj, 1])
    out = {f: np.array(cg[("boxlib", f)])[:, :, 0].T for f in fields}
    return (i0 + 0.5 + np.arange(ni)) * dx, (j0 + 0.5 + np.arange(nj)) * dx, out, dx


def one(case):
    try:
        fi = os.path.join(PELE, case, "input.case")
        inp = open(fi if os.path.exists(fi) else os.path.join(PELE, case, "input.c1")).read()
        g = lambda k: re.search(r"^%s\s*=\s*(\S+)" % re.escape(k), inp, re.M).group(1)
        P = float(g("prob.P_mean")); p = 1 if P < 2e5 else int(round(P / 1e5)); Tu = float(g("prob.T_in"))
        poly = sL_poly(*TABLE[(p, int(round(Tu)))])
        pf = sorted(glob.glob(os.path.join(PELE, case, "plt?????")))[-1]
        ds = yt.load(pf)
        W, H = float(ds.domain_width[0]), float(ds.domain_width[1])
        dL = 317e-6 * W / 0.075
        x, y, f, _ = cgrid(ds, 1, (0.0, 0.0), (W, H), ["temp"])
        T = f["temp"]
        j, i = np.unravel_index(np.argmin(np.where(T > 1000.0, y[:, None] + 0 * x[None, :], 1e9)), T.shape)
        lev = ds.max_level
        x, y, f, dx = cgrid(ds, lev, (x[i] - 8 * dL, y[j] - 6 * dL), (x[i] + 8 * dL, y[j] + 10 * dL), F)
        T, u, v, rho = f["temp"], f["x_velocity"], f["y_velocity"], f["density"]
        pts = np.concatenate([l for l in contourpy.contour_generator(x, y, T).lines(1000.0) if len(l) > 10])
        x0, y0 = pts[np.argmin(pts[:, 1])]
        # fields for the stretch (arrays indexed [j, i] = [y, x])
        Ty, Tx = np.gradient(T, dx)
        gm = np.hypot(Tx, Ty) + 1e-30
        nx, ny = -Tx / gm, -Ty / gm
        kappa = np.gradient(nx, dx, axis=1) + np.gradient(ny, dx, axis=0)
        uy, ux = np.gradient(u, dx); vy, vx = np.gradient(v, dx)
        div = ux + vy
        nn = nx * nx * ux + nx * ny * uy + ny * nx * vx + ny * ny * vy
        at = div - nn
        sd = -(u * nx + v * ny)
        it = lambda a, xx, yy: float(RegularGridInterpolator((y, x), a)((yy, xx)))
        ic = int(np.argmin(np.abs(x - x0)))
        col = lambda a: a[:, ic]
        # vertical through the tip: unburnt edge, consumption speed
        Tc, yc = col(T), y
        ju = int(np.where(Tc[: np.searchsorted(yc, y0)] < Tu + 10.0)[0][-1])     # last unburnt cell below the tip
        yu = float(np.interp(Tu + 10.0, Tc[ju:ju + 2], yc[ju:ju + 2]))
        wO2 = -col(f["I_R(O2)"])
        kp = int(np.argmax(wO2))
        ke = kp + 1
        while ke < len(wO2) - 1 and not (wO2[ke] <= wO2[ke - 1] and wO2[ke] <= wO2[ke + 1]):
            ke += 1
        j0 = max(ju - int(2 * dL / dx), 0)
        rho_u, YO2u, zu = col(rho)[ju], col(f["Y(O2)"])[ju], col(f["mixture_fraction"])[ju]
        sc = float(np.trapezoid(wO2[j0:ke + 1], yc[j0:ke + 1]) / (rho_u * YO2u))
        o2_left = float(col(f["Y(O2)"])[ke] / YO2u)
        sLz = float(np.polyval(poly, zu))
        r = dict(case=case, plt=os.path.basename(pf), p=p, Tu=Tu, dL=dL, dx=dx, x_tip=float(x0), y_tip=float(y0),
                 z_u=float(zu), phi_u=float(zu * (1 - 0.02830074092068505) / (0.02830074092068505 * (1 - zu))),
                 sL=sLz, s_c=sc, sc_over_sL=sc / sLz, O2_left_end=o2_left,
                 window_dL=float((yc[ke] - yc[j0]) / dL))
        for tag, yy in (("u", yu), ("1000", float(y0))):
            K_s, K_k = it(at, x0, yy), it(sd, x0, yy) * it(kappa, x0, yy)
            r.update({"kappa_dL_" + tag: it(kappa, x0, yy) * dL, "sd_" + tag: it(sd, x0, yy),
                      "sd_rho_" + tag: it(sd * rho, x0, yy) / rho_u,
                      "Ka_s_" + tag: K_s * dL / sLz, "Ka_k_" + tag: K_k * dL / sLz, "Ka_" + tag: (K_s + K_k) * dL / sLz})
        # (2) eps* seen by the tip: horizontal lines below the tip, finest level (fextract, as eps_measured.py; a long
        # 2D covering grid of a few rows fails in yt)
        for k in (2, 4, 8):
            yl = y0 - k * dL
            out = os.path.join(PELE, case, "slice_up%d.dat" % k)
            subprocess.run([FX, "-d", "0", "-y", "%.8e" % yl, "-v", "mixture_fraction", "-s", out, pf], check=True,
                           stdout=subprocess.DEVNULL)
            xl, zl = np.loadtxt(out, unpack=True)
            xl, iu = np.unique(xl, return_index=True); zl = zl[iu]
            xd = xl / dL
            xf = np.arange(xd[0], xd[-1], 0.002)
            zf = PchipInterpolator(xd, zl)(xf)
            ok = (zf > 0.005) & (zf < phi2z(2.0))
            r["eps_up%d" % k] = float(eps_fit(xf[ok], np.polyval(poly, zf[ok]))[0])
            r["dplus_up%d" % k] = float(1.0 / np.max(np.gradient(zl, xd)))
            r["z_up%d" % k] = float(np.interp(x0 / dL, xd, zl))
        # reference flame speed: z on the vertical through the tip 8 deltaL upstream, before the preferential diffusion
        # of H2 enriches the preheat zone (z_u at T = Tu + 10 K is far richer than the layer)
        zr = r["z_up8"]
        r["sL8"] = float(np.polyval(poly, zr))
        r["phi8"] = float(zr * (1 - 0.02830074092068505) / (0.02830074092068505 * (1 - zr)))
        r["sc_over_sL8"] = sc / r["sL8"]
        for tag in ("u", "1000"):
            for q in ("Ka_", "Ka_s_", "Ka_k_"):
                r[q + tag + "_8"] = r[q + tag] * sLz / r["sL8"]
        return r
    except Exception as e:                                       # noqa: BLE001
        return dict(case=case, error=repr(e))


if __name__ == "__main__":
    with Pool(int(sys.argv[2])) as pool:
        res = pool.map(one, sys.argv[3:])
    for r in res:
        if "error" in r:
            print(r["case"], "ERROR", r["error"]); continue
        print("%-15s phi8=%.2f sc/sL8=%.3f phi_u=%.2f sc/sL=%.3f O2left=%.2f Ka_u=%.3f (s %.3f k %.3f) Ka_1000=%.3f kdL=%.3f "
              "eps_up2/4/8=%.3f/%.3f/%.3f" % (r["case"], r["phi8"], r["sc_over_sL8"], r["phi_u"], r["sc_over_sL"], r["O2_left_end"], r["Ka_u"],
                                              r["Ka_s_u"], r["Ka_k_u"], r["Ka_1000"], r["kappa_dL_1000"],
                                              r["eps_up2"], r["eps_up4"], r["eps_up8"]), flush=True)
    json.dump(res, open(sys.argv[1], "w"), indent=1)
