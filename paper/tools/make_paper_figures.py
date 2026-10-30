"""Figures of the JFM Rapids manuscript (paper/figures/*.pdf), compact multi-panel layout.
Fig. 1: (a) heat-release field of a steady triple flame (PeleLMeX, case C4, 1 atm, mixture-averaged transport)
        with mixture-fraction isolines and streamlines;
        (b) normalised laminar flame speed sL/s_m versus phi at the six conditions (Cantera, R0 table);
        (c) s(x)/s_m across erf and tanh layers of equal Delta, with the osculating parabolas (tip radius R).
Fig. 2: data of Chaplet et al. (2026) versus 1/Delta+ (a) and eps* (b), with the C3-C5 extrapolation test.
Fig. 3: sT/sTth versus Delta+ (a) and eps* (b) for the two transport formulations ('ma_*': mixture-averaged,
        'sc_*': constant Schmidt numbers, see data/results_current.json) at all conditions, with the fitted lines;
        (c) slope ratio b/b(1 atm) against lambda Ma_c and lambda Ma_u; (d) measured against predicted slope.
        The fits are written to data/fits.json.
Data: data/results_current.json, theory/sL_cantera_phi.csv, $RUNS/C4 (last plotfile).
Usage: python tools/make_paper_figures.py"""
import csv, glob, json, os
import numpy as np
from scipy.signal import fftconvolve
from scipy.special import erf
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter, FixedLocator, FixedFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
ROOT = os.path.dirname(PAPER)
OUT = os.path.join(PAPER, "figures")
os.makedirs(OUT, exist_ok=True)
D = json.load(open(os.path.join(PAPER, "data", "results_current.json")))
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6, "lines.linewidth": 1.0, "savefig.bbox": "tight",
                     "font.family": "serif", "mathtext.fontset": "cm", "xtick.direction": "in",
                     "ytick.direction": "in", "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "legend.handlelength": 1.5})
W2 = 5.3                      # JFM text width [in]
zst = 0.02830074092068505
phi2z = lambda phi: phi * zst / (1.0 - zst + phi * zst)
dL1 = 317e-6


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=170)
    plt.close(fig)


def lab(ax, s):
    ax.text(0.02, 0.97, s, transform=ax.transAxes, va="top", ha="left", fontsize=8)


def r0_table(cond, dx):
    rows = [r for r in csv.DictReader(open(os.path.join(ROOT, "theory", "sL_cantera_phi.csv"))) if r["condition"] == cond]
    d = {}
    for r in rows:
        d[float(r["phi"])] = float(r["sL_cantera_m_s"])
    phi = np.array(sorted(d))
    return phi, np.array([d[p] for p in phi])


def sL_poly(cond, dx):
    phi, sl = r0_table(cond, dx)
    return np.polyfit(phi2z(phi), sl, 4)


def layer(x, d0p, dp):
    """tanh inlet profile of thickness d0p diffused by a Gaussian to the thickness dp (deltaL units)."""
    h = x[1] - x[0]
    z = 0.5 * (1 + np.tanh(2 * x / d0p))
    dd = np.sqrt(max(dp ** 2 - d0p ** 2, 0.0))
    if dd > 0:
        sd = dd / np.sqrt(2 * np.pi)
        k = np.exp(-0.5 * (np.arange(-int(8 * sd / h), int(8 * sd / h) + 1) * h / sd) ** 2)
        z = fftconvolve(z, k / k.sum(), mode="same")
    return z


def eps_star(P, d0p, dp):
    x = np.linspace(-80, 80, 160001)
    h = x[1] - x[0]
    z = layer(x, d0p, dp)
    s = np.polyval(P, z)
    ok = (z > 0.005) & (z < phi2z(2.0))
    m = np.argmax(np.where(ok, s, -1))
    return np.sqrt(abs(np.gradient(np.gradient(s, h), h)[m]) / s[m])


P1 = sL_poly("p1_T300", 25.0)
P10 = sL_poly("p10_T579", 2.0)
P20 = sL_poly("p20_T706", 0.92)

# =========================================================== Fig. 1
fig = plt.figure(figsize=(W2, 2.35))
gs = fig.add_gridspec(1, 3, width_ratios=[1.08, 1, 1], wspace=0.48)
a = fig.add_subplot(gs[0])
import yt
from matplotlib import patheffects as pe
yt.set_log_level(50)
pf = sorted(glob.glob(os.path.join(os.environ.get("RUNS", os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "cases")), "C4", "plt?????")))[-1]
ds = yt.load(pf)
lev = 3
dxl = float(ds.domain_width[0]) / (ds.domain_dimensions[0] * 2 ** lev)
x0, x1, y0w, y1w = 0.0375 - 14 * dL1, 0.0375 + 14 * dL1, 0.0187 - 9 * dL1, 0.0187 + 31 * dL1
i0, j0 = int(x0 / dxl), int(y0w / dxl)
ni, nj = int((x1 - x0) / dxl), int((y1w - y0w) / dxl)
cg = ds.covering_grid(lev, left_edge=[i0 * dxl, j0 * dxl, 0.0], dims=[ni, nj, 1])
g = lambda f: np.array(cg[("boxlib", f)])[:, :, 0].T
hr, z, u, v, T = g("HeatRelease"), g("mixture_fraction"), g("x_velocity"), g("y_velocity"), g("temp")
xs = (i0 + 0.5 + np.arange(ni)) * dxl; ys = (j0 + 0.5 + np.arange(nj)) * dxl
sh = T.shape
jt, it = np.unravel_index(np.argmin(np.where(T > 1000, ys[:, None] + 0 * xs[None, :], 1.0)), sh)
X, Y = (xs - xs[it]) / dL1, (ys - ys[jt]) / dL1
sx, sy = (abs(X) < 14), (Y > -8) & (Y < 30)
Xs, Ys = X[sx], Y[sy]
hrn = hr[np.ix_(sy, sx)] / hr.max()
pc = a.pcolormesh(Xs, Ys, np.clip(hrn, 1e-3, 1), norm=LogNorm(1e-3, 1), cmap="inferno", shading="auto", rasterized=True)
a.streamplot(Xs, Ys, u[np.ix_(sy, sx)], v[np.ix_(sy, sx)], color="#9aa5b1", linewidth=0.3, density=0.6, arrowsize=0.35)
zz_ = z[np.ix_(sy, sx)]
cl = a.contour(Xs, Ys, zz_, levels=[phi2z(1.0), phi2z(1.6)], colors=["w", "#5ec8ff"], linewidths=0.8)
a.contour(Xs, Ys, zz_, levels=[phi2z(0.5), phi2z(3.0)], colors="w", linewidths=0.45, linestyles="--")
halo = [pe.withStroke(linewidth=1.6, foreground="k")]
for txt, xx_, yy_, col in ((r"$\phi=1$", -6.3, 21.0, "w"), (r"$\phi=1.6$", 0.6, 21.0, "#5ec8ff")):
    a.text(xx_, yy_, txt, color=col, fontsize=6, ha="left", va="center", path_effects=halo)
a.plot([0], [0], "o", ms=2.5, mfc="#5ec8ff", mec="k", mew=0.4)
a.set_xlim(-14, 14); a.set_ylim(-8, 30); a.set_aspect("equal")
a.set_xlabel(r"$(x-x_{tip})/\delta_L$"); a.set_ylabel(r"$(y-y_{tip})/\delta_L$")
a.text(0.04, 0.97, r"($a$)", transform=a.transAxes, va="top", color="w", fontsize=8, path_effects=halo)
cax = a.inset_axes([0.0, 1.04, 1.0, 0.045])
cb = fig.colorbar(pc, cax=cax, orientation="horizontal"); cb.set_ticks([1e-3, 1e-2, 1e-1, 1])
cb.ax.xaxis.set_ticks_position("top"); cb.ax.xaxis.set_label_position("top")
cb.ax.tick_params(labelsize=5.5, length=2, width=0.5, pad=1); cb.outline.set_linewidth(0.5)
cb.set_label(r"heat release $/$ max", fontsize=6, labelpad=2)
b = fig.add_subplot(gs[1])
blues = ["#0b3c5d", "#1f6fa8", "#4fa3d8", "#9cc9ea"]
conds = (("p1_T300", 25.0, "1 bar, 300 K", blues[0], "-"), ("p5_T475", 4.1, "5 bar, 475 K", blues[1], "-"),
         ("p10_T579", 2.0, "10 bar, 579 K", blues[2], "-"), ("p20_T706", 0.92, "20 bar, 706 K", blues[3], "-"),
         ("p10_T300", 1.6, "10 bar, 300 K", "#d1495b", "--"), ("p1_T579", 34.0, "1 bar, 579 K", "#e8a33d", "--"))
for cnd, dx, lb, col, ls_ in conds:
    phi, sl = r0_table(cnd, dx)
    pp = np.polyfit(phi, sl, 4); ph = np.linspace(1.0, 2.0, 100); sm = np.polyval(pp, ph).max()
    b.plot(ph, np.polyval(pp, ph) / sm, ls_, color=col, lw=0.9, label=lb)
    b.plot(phi, sl / sm, "o", ms=2.0, color=col, mfc="w", mew=0.6)
b.set_xlabel(r"$\phi$"); b.set_ylabel(r"$s_L/s_m$"); b.set_ylim(0.6, 1.04); b.set_xlim(0.95, 2.05)
b.legend(frameon=False, fontsize=5.6, loc="lower right", handlelength=1.6, labelspacing=0.25)
lab(b, r"($b$)")
c = fig.add_subplot(gs[2])
xx = np.linspace(-3, 3, 60001)
for prof, zz, ls_, ha, dxl_ in (("erf", 0.5 * (1 + erf(np.sqrt(np.pi) * xx)), "-", "left", 0.012),
                                ("tanh", 0.5 * (1 + np.tanh(2 * xx)), "--", "right", -0.012)):
    s = np.polyval(P1, zz); ok = (zz > phi2z(1.0)) & (zz < phi2z(2.0)); mm = np.argmax(np.where(ok, s, -1))
    s2 = np.gradient(np.gradient(s, xx[1] - xx[0]), xx[1] - xx[0])[mm]
    c.plot(xx[ok], s[ok] / s[mm], "k" + ls_, lw=0.9)
    xp = np.linspace(xx[mm] - 0.07, xx[mm] + 0.07, 30)
    c.plot(xp, 1 + 0.5 * s2 / s[mm] * (xp - xx[mm]) ** 2, ":", color="#d1495b", lw=1.3)
    c.text(xx[mm] + dxl_, 1.018, r"%s" % prof + "\n" + r"$C=%.1f$" % np.sqrt(abs(s2) / s[mm]), ha=ha, va="bottom",
           fontsize=6.2, linespacing=1.1)
c.set_xlim(-1.0, -0.55); c.set_ylim(0.74, 1.09)
c.set_xlabel(r"$x/\Delta$"); c.set_ylabel(r"$s_L(z(x))/s_m$")
c.text(-0.605, 0.835, "osculating\nparabola", color="#d1495b", fontsize=5.6, ha="center", va="top", linespacing=1.0)
lab(c, r"($c$)")
save(fig, "fig1")

# =========================================================== Fig. 2
ch = D["chaplet"]; names = list(ch)
d0 = np.array([ch[n][0] for n in names]); dp = np.array([ch[n][1] for n in names]); r = np.array([ch[n][2] for n in names])
eps = np.array([eps_star(P1, a_, b_) for a_, b_ in zip(d0, dp)])
json.dump({n: dict(d0p=a_, dplus=b_, ratio=c_, eps=e_) for n, a_, b_, c_, e_ in zip(names, d0, dp, r, eps)},
          open(os.path.join(PAPER, "data", "chaplet_eps.json"), "w"), indent=1)
fig, ax = plt.subplots(1, 2, figsize=(W2, 2.0), sharey=True, gridspec_kw=dict(wspace=0.08))
fitC = {}
CB, CR, CO = "#1f6fa8", "#d1495b", "#e8a33d"
OFF = {"C1": (0, -10, "center"), "C2": (3, 4, "left"), "C3": (-4, -9, "right"), "C4": (-4, -9, "right"),
       "C5": (-4, -9, "right"), "C6": (-5, 4, "right")}
# (a) functions of Delta+ fitted on C3-C5 (two parameters each), Delta+ on a log axis
a = ax[0]
dg = np.exp(np.linspace(np.log(1.8), np.log(40), 200))
for f, lbl, col, ls_ in ((lambda x: 1 / x, r"$a+b/\Delta^+$", CB, "--"), (lambda x: x ** -0.5, r"$a+b\,(\Delta^+)^{-1/2}$", CO, "-.")):
    bb_, aa_ = np.polyfit(f(dp[2:5]), r[2:5], 1)
    a.plot(dg, aa_ + bb_ * f(dg), ls_, color=col, lw=0.9, label=lbl + "  (C6: %+.3f)" % (aa_ + bb_ * f(dp[5]) - r[5]))
    a.plot(dp[5], aa_ + bb_ * f(dp[5]), "x", color=col, ms=5, mew=1.0, zorder=3)
a.plot(dp, r, "o", color="#1a1a1a", ms=3.6, mfc="w", mew=0.8, zorder=4)
for n_, x_, y_ in zip(names, dp, r):
    dx_, dy_, ha = {"C1": (-4, 4, "right"), "C2": (0, 5, "center"), "C3": (4, -8, "left"), "C4": (4, -8, "left"),
                    "C5": (4, -8, "left"), "C6": (4, 3, "left")}[n_]
    a.annotate(n_, (x_, y_), xytext=(dx_, dy_), textcoords="offset points", fontsize=6, ha=ha)
a.set_xscale("log"); a.set_xlim(1.8, 40)
a.xaxis.set_major_locator(FixedLocator([2, 3, 5, 10, 20, 40])); a.xaxis.set_major_formatter(FixedFormatter(["2", "3", "5", "10", "20", "40"]))
a.xaxis.set_minor_formatter(NullFormatter())
a.set_xlabel(r"$\Delta^+$"); a.legend(frameon=False, fontsize=5.6, loc="lower right", handlelength=1.8,
                                       title="fits on C3–C5 (prediction error on C6)", title_fontsize=5.6)
a.text(0.03, 0.96, r"($a$)", transform=a.transAxes, va="top", ha="left", fontsize=8)
# (b) eps*
a = ax[1]; xv = eps
b3, c3 = np.polyfit(xv[2:5], r[2:5], 1); b4, c4 = np.polyfit(xv[2:6], r[2:6], 1)
xg = np.linspace(0, xv.max() * 1.08, 50)
a.plot(xg, c3 + b3 * xg, "--", color=CB, lw=0.9, label="fit on C3–C5  (C6: %+.3f)" % (c3 + b3 * xv[5] - r[5]))
a.plot(xg, c4 + b4 * xg, "-", color=CR, lw=0.9, label="fit on C3–C6")
a.plot(xv[5], c3 + b3 * xv[5], "x", color=CB, ms=5, mew=1.0, zorder=3)
a.plot(xv, r, "o", color="#1a1a1a", ms=3.6, mfc="w", mew=0.8, zorder=4)
for n_, x_, y_ in zip(names, xv, r):
    dx_, dy_, ha = OFF[n_]
    a.annotate(n_, (x_, y_), xytext=(dx_, dy_), textcoords="offset points", fontsize=6, ha=ha)
a.set_xlim(0, xv.max() * 1.08); a.set_xlabel(r"$\varepsilon^*=\delta_L/R$")
a.legend(frameon=False, fontsize=5.6, loc="upper right", handlelength=1.8)
a.text(0.03, 0.04, r"($b$)", transform=a.transAxes, va="bottom", ha="left", fontsize=8)
fitC["b"] = (c4, b4)
ax[0].set_ylabel(r"$s_T/s_{T,th}$"); ax[0].set_ylim(0.52, 1.12)
save(fig, "fig2")
A0, B0 = fitC["b"]

# =========================================================== Fig. 3
th = D["sTth"]
P20 = sL_poly("p20_T706", 0.92)
POLY = {"1bar": P1, "10bar": P10, "20bar": P20, "1bar579K": sL_poly("p1_T579", 34.0), "10bar300K": sL_poly("p10_T300", 1.6),
        "5bar": sL_poly("p5_T475", 4.1)}
NAME = {"sc": "constant Sc", "ma": "mixture-averaged"}
SETS = [("sc", p, NAME["sc"], "o") for p in ("1bar", "10bar", "20bar", "10bar300K", "1bar579K")] + \
       [("ma", p, NAME["ma"], "^") for p in ("1bar", "5bar", "10bar", "20bar", "10bar300K", "1bar579K")]
COL = {"1bar": "#1a1a1a", "5bar": "#0b3d63", "10bar": "#1f6fa8", "20bar": "#7cc0ea", "1bar579K": "#e8a33d", "10bar300K": "#d1495b"}
EDGE = dict(COL, **{"20bar": "#3b8fc9"})
PLAB = {"1bar": "1 atm", "5bar": "5 bar", "10bar": "10 bar", "20bar": "20 bar", "1bar579K": "1 bar, 579 K", "10bar300K": "10 bar, 300 K"}
data, fits = {}, {}
for code, p, name, mk in SETS:
    rows = []
    for n, v_ in D["%s_%s" % (code, p)].items():
        if v_.get("dplus"):
            rows.append((n, v_["dplus"], eps_star(POLY[p], v_["d0p"], v_["dplus"]), v_["sT"] / th["%s_%s" % (code, p)],
                         bool(v_.get("prov"))))
    data[(code, p)] = rows
    e_ = np.array([r_[2] for r_ in rows]); s_ = np.array([r_[3] for r_ in rows])
    (bb, aa), cov = np.polyfit(e_, s_, 1, cov=True)
    fits["%s_%s" % (code, p)] = dict(A=float(aa), b=float(-bb), sb=float(np.sqrt(cov[0, 0])), n=len(e_),
                                     rms=float(np.sqrt(np.mean((s_ - aa - bb * e_) ** 2))), eps_max=float(e_.max()))
for code in ("sc", "ma"):
    b1 = fits["%s_1bar" % code]
    for p in ("5bar", "10bar", "20bar", "10bar300K", "1bar579K"):
        k = "%s_%s" % (code, p)
        if k in fits:
            r_ = fits[k]["b"] / b1["b"]
            fits[k]["ratio"] = r_
            fits[k]["sratio"] = r_ * np.hypot(fits[k]["sb"] / fits[k]["b"], b1["sb"] / b1["b"])
json.dump(fits, open(os.path.join(PAPER, "data", "fits.json"), "w"), indent=1)

from matplotlib.lines import Line2D
fig, ax = plt.subplots(2, 2, figsize=(W2, 3.95), gridspec_kw=dict(wspace=0.3, hspace=0.45))
a, b, c, d = ax.ravel()
a.plot(dp, r, "s", color="#6b7785", ms=3, mfc="w", mew=0.6)
b.plot(eps, r, "s", color="#6b7785", ms=3, mfc="w", mew=0.6)
for code, p, name, mk in SETS:
    for n, dpp, ee, rr, prov in data[(code, p)]:
        kw = dict(ms=3.7, color=EDGE[p], mfc="w" if (prov or code == "ma") else COL[p], mew=0.75)
        a.plot(dpp, rr, mk, **kw)
        b.plot(ee, rr, mk, **kw)
    f_ = fits["%s_%s" % (code, p)]
    eg = np.linspace(0.1, f_["eps_max"] * 1.05, 20)
    b.plot(eg, f_["A"] - f_["b"] * eg, "-" if code == "sc" else "--", color=EDGE[p], lw=0.75)
H = [Line2D([], [], ls="", marker="o", ms=4, color=EDGE[p], mfc=COL[p], label=PLAB[p])
     for p in ("1bar", "5bar", "10bar", "20bar", "10bar300K", "1bar579K")]
H += [Line2D([], [], ls="-", marker="o", ms=4, color="0.35", mfc="0.35", lw=0.8, label="constant Sc"),
      Line2D([], [], ls="--", marker="^", ms=4, color="0.35", mfc="w", lw=0.8, label="mixture-averaged"),
      Line2D([], [], ls="", marker="s", ms=3.5, color="#6b7785", mfc="w", label="Chaplet et al. (2026)")]
fig.canvas.draw()
fig.legend(handles=H, frameon=False, fontsize=5.8, loc="lower center", ncol=5, columnspacing=0.9, handletextpad=0.3,
           bbox_to_anchor=(0.5, a.get_position().y1 + 0.005), handlelength=1.8)
a.set_xscale("log"); a.set_xlim(1.6, 40); a.set_ylim(0.62, 1.12)
a.xaxis.set_major_locator(FixedLocator([2, 3, 5, 10, 20, 40])); a.xaxis.set_major_formatter(FixedFormatter(["2", "3", "5", "10", "20", "40"]))
a.xaxis.set_minor_formatter(NullFormatter())
a.set_xlabel(r"$\Delta^+$"); a.set_ylabel(r"$s_T/s_{T,th}$"); lab(a, r"($a$)")
b.set_xlim(0, 4.2); b.set_ylim(0.62, 1.12); b.set_xlabel(r"$\varepsilon^*$"); b.set_ylabel(r"$s_T/s_{T,th}$")
b.text(0.03, 0.04, r"($b$)", transform=b.transAxes, va="bottom", ha="left", fontsize=8)
pr = D["prediction_b_ratio"]
# (c) isentropic path against pressure (continuous prediction, data/markstein_sweep.json) | single-effect conditions
SW = sorted((r_["p"], r_["b_c"]) for r_ in json.load(open(os.path.join(PAPER, "data", "markstein_sweep.json")))
            if r_["mech"] == "sandiego" and "error" not in r_)
sw_p = np.array([s_[0] for s_ in SW]); sw_r = np.array([s_[1] for s_ in SW]) / SW[0][1]
pos = c.get_position()
c.set_position([pos.x0, pos.y0, pos.width * 0.60, pos.height])
c2 = fig.add_axes([pos.x0 + pos.width * 0.66, pos.y0, pos.width * 0.34, pos.height], sharey=c)
pg = np.linspace(1, 20, 200)
c.plot(pg, np.interp(pg, sw_p, sw_r), "-", color="#1a1a1a", lw=1.2, label=r"$\lambda\,Ma_c$")
for p, pp in (("5bar", 5.0), ("10bar", 10.0), ("20bar", 20.0)):
    c.plot([pp - 0.9, pp + 0.9], [pr[p]["Ma_u"]] * 2, "--", color="0.6", lw=1.1, label=r"$\lambda\,Ma_u$" if p == "10bar" else None)
PISO = {"5bar": 5.0, "10bar": 10.0, "20bar": 20.0}
for p, pp in PISO.items():
    for code, mk, dx in (("sc", "o", -0.45), ("ma", "^", 0.45)):
        kk = "%s_%s" % (code, p)
        if kk in fits:
            c.errorbar(pp + dx, fits[kk]["ratio"], yerr=fits[kk]["sratio"], fmt=mk, color=EDGE[p], ms=3.8,
                       mfc="w" if code == "ma" else COL[p], mew=0.75, elinewidth=0.7, capsize=1.6)
c.plot(1, 1, "o", color="#1a1a1a", ms=3.2)
c.axhline(1, color="0.6", lw=0.5, ls=":")
c.set_xlim(0, 21.5); c.set_ylim(0.5, 2.45); c.set_xticks([1, 5, 10, 15, 20])
c.set_xlabel(r"$p$ (bar), isentropic path", fontsize=7)
c.set_ylabel(r"$b/b(1\,\mathrm{atm})$"); lab(c, r"($c$)")
c.legend(frameon=False, fontsize=6, loc="lower right", ncol=1, handlelength=1.8, borderaxespad=0.3)
CATS = [("10bar300K", "10 bar\n300 K"), ("1bar579K", "1 bar\n579 K")]
for k_, (p, lb) in enumerate(CATS):
    c2.plot([k_ - 0.32, k_ + 0.32], [pr[p]["Ma_c"]] * 2, "-", color="#1a1a1a", lw=1.3)
    c2.plot([k_ - 0.32, k_ + 0.32], [pr[p]["Ma_u"]] * 2, "--", color="0.6", lw=1.1)
    for code, mk, dx in (("sc", "o", -0.12), ("ma", "^", 0.12)):
        kk = "%s_%s" % (code, p)
        if kk in fits:
            c2.errorbar(k_ + dx, fits[kk]["ratio"], yerr=fits[kk]["sratio"], fmt=mk, color=EDGE[p], ms=3.8,
                        mfc="w" if code == "ma" else COL[p], mew=0.75, elinewidth=0.7, capsize=1.6)
c2.axhline(1, color="0.6", lw=0.5, ls=":")
c2.set_xlim(-0.55, len(CATS) - 0.45); c2.set_xticks(range(len(CATS)))
c2.set_xticklabels([lb for _, lb in CATS], fontsize=5.8)
c2.tick_params(labelleft=False)
c2.text(0.5, 1.02, "single effect", transform=c2.transAxes, ha="center", va="bottom", fontsize=6, color="0.35")
PB = D["prediction_b"]
bg = np.linspace(0, 0.3, 10)
ratio = np.mean([fits[k]["b"] / PB[k] for k in fits])
# 30/10/2026 (review point 3): mean-ratio line and +-15 % band removed, equality only
d.plot(bg, bg, ":", color="0.4", lw=0.7)
for code, p, name, mk in SETS:
    k = "%s_%s" % (code, p)
    d.plot(PB[k], fits[k]["b"], mk, color=EDGE[p], ms=4.2, mfc="w" if code == "ma" else COL[p], mew=0.8, zorder=3)
d.text(0.205, 0.183, "1:1", fontsize=6, color="0.35", rotation=36)
d.set_xlim(0, 0.27); d.set_ylim(0, 0.29)
d.set_xlabel(r"$\lambda\,Ma_{c,\delta}$ (prediction)"); d.set_ylabel(r"$b$ (simulations)")
d.text(0.97, 0.04, r"($d$)", transform=d.transAxes, va="bottom", ha="right", fontsize=8)
print("parity mean ratio b/(lambda Ma_c) = %.3f" % ratio, {k: round(fits[k]["b"] / PB[k], 3) for k in fits})
save(fig, "fig3")
print("eps* Chaplet:", dict(zip(names, np.round(eps, 3))))
for k, v_ in fits.items():
    print("%-12s" % k, {kk: round(vv, 3) if isinstance(vv, float) else vv for kk, vv in v_.items()})
