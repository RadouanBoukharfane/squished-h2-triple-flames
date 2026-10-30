"""Illustrative figures of the PeleLMeX runs (talks, highlights; not the manuscript), paper/highlights/*.png:
H1 pressure : same Delta+ ~ 5 at 1, 10 and 20 bar (heat release in deltaL units, 1000 K isotherm, stoichiometric line)
H2 squishing: C1, C4, C6 at 1 atm
H3 transport: the same inlet layer, cold, mixture-averaged against constant Schmidt numbers (z field and Delta(y))
H4 inlet    : tanh layer held 3 mm from the inlet against an erf layer held far from it
Usage (compute node): python tools/highlights.py"""
import glob, json, os
import numpy as np
import yt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

yt.set_log_level(50)
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PELE = os.environ.get("RUNS", os.path.join(ROOT, "cases"))
OUT = os.path.join(os.path.dirname(HERE), "highlights")
os.makedirs(OUT, exist_ok=True)
ZST = 0.02830074092068505
plt.rcParams.update({"font.size": 11, "font.family": "serif", "mathtext.fontset": "cm", "axes.linewidth": 0.7,
                     "savefig.bbox": "tight", "savefig.dpi": 220})


def field(case, fields, xw=(-12, 12), yw=(-6, 24), lev=None, cold=False):
    pf = sorted(glob.glob(os.path.join(PELE, case, "plt?????")))[-1]
    ds = yt.load(pf)
    W, H = float(ds.domain_width[0]), float(ds.domain_width[1])
    dL = 317e-6 * W / 0.075
    lev = ds.max_level if lev is None else min(lev, ds.max_level)
    dx = W / (ds.domain_dimensions[0] * 2 ** lev)
    if cold:
        xt, yt_ = 0.5 * W, 0.0
    else:   # leading edge from a level-1 pass
        dx1 = W / (ds.domain_dimensions[0] * 2)
        cg = ds.covering_grid(1, left_edge=[0, 0, 0], dims=[ds.domain_dimensions[0] * 2, ds.domain_dimensions[1] * 2, 1])
        T = np.array(cg[("boxlib", "temp")])[:, :, 0].T
        ys = (np.arange(T.shape[0]) + 0.5) * dx1; xs = (np.arange(T.shape[1]) + 0.5) * dx1
        j, i = np.unravel_index(np.argmin(np.where(T > 1000, ys[:, None] + 0 * xs[None, :], 1e9)), T.shape)
        xt, yt_ = xs[i], ys[j]
    x0, x1 = max(xt + xw[0] * dL, 0), min(xt + xw[1] * dL, W)
    y0, y1 = max(yt_ + yw[0] * dL, 0), min(yt_ + yw[1] * dL, H)
    i0, j0 = int(x0 / dx), int(y0 / dx)
    ni, nj = int((x1 - x0) / dx), int((y1 - y0) / dx)
    cg = ds.covering_grid(lev, left_edge=[i0 * dx, j0 * dx, 0], dims=[ni, nj, 1])
    F = {f: np.array(cg[("boxlib", f)])[:, :, 0].T for f in fields}
    X = ((i0 + 0.5 + np.arange(ni)) * dx - xt) / dL
    Y = ((j0 + 0.5 + np.arange(nj)) * dx - yt_) / dL
    return X, Y, F, dict(dL=dL, y_tip=yt_, xt=xt)


def flame_panel(ax, case, title, sub, xw=(-12, 12), yw=(-6, 24), show_inlet=False):
    X, Y, F, m = field(case, ["HeatRelease", "temp", "mixture_fraction", "x_velocity", "y_velocity"], xw, yw)
    hr = F["HeatRelease"]
    pc = ax.pcolormesh(X, Y, np.clip(hr / hr.max(), 1e-3, 1), norm=LogNorm(1e-3, 1), cmap="inferno", shading="auto",
                       rasterized=True)
    ax.contour(X, Y, F["temp"], levels=[1000], colors="#7fd3ff", linewidths=0.9)
    ax.contour(X, Y, F["mixture_fraction"], levels=[ZST], colors="w", linewidths=0.6, linestyles="--")
    ax.streamplot(X, Y, F["x_velocity"], F["y_velocity"], color="0.7", linewidth=0.4, density=0.8, arrowsize=0.5)
    if show_inlet:
        yin = -m["y_tip"] / m["dL"]
        ax.axhspan(yin - 3, yin, color="#3a3a3a"); ax.axhline(yin, color="#f4a261", lw=2.5)
        ax.text(X[-1] - 0.5, yin + 0.6, "inlet", color="#f4a261", fontsize=11, ha="right")
    ax.set_xlim(X[0], X[-1]); ax.set_ylim(min(Y[0], yw[0]) if show_inlet else Y[0], Y[-1]); ax.set_aspect("equal")
    ax.set_title(title, fontsize=12)
    ax.text(0.03, 0.97, sub, transform=ax.transAxes, va="top", ha="left", color="w", fontsize=10,
            bbox=dict(facecolor="k", alpha=0.45, lw=0, pad=3))
    ax.set_xlabel(r"$(x-x_{tip})/\delta_L$")
    return pc


def three(cases, name, suptitle, **kw):
    fig, ax = plt.subplots(1, len(cases), figsize=(4.0 * len(cases), 5.2), sharey=True)
    for a, (case, title, sub) in zip(ax, cases):
        pc = flame_panel(a, case, title, sub, **kw)
    ax[0].set_ylabel(r"$(y-y_{tip})/\delta_L$")
    cb = fig.colorbar(pc, ax=ax, shrink=0.8, pad=0.015); cb.set_label("heat release / max")
    fig.canvas.draw()
    fig.suptitle(suptitle, fontsize=13, y=ax[0].get_position().y1 + 0.09)
    fig.savefig(os.path.join(OUT, name + ".png")); plt.close(fig)
    print("wrote", name, flush=True)


# H1 pressure
three([("C4", "1 atm, 300 K", r"$\Delta^+=4.9$" "\n" r"$s_T/s_{T,th}=0.93$"),
       ("B1-p10-D5b", "10 bar, 579 K", r"$\Delta^+=5.1$" "\n" r"$s_T/s_{T,th}=0.87$"),
       ("B1-p20-D5", "20 bar, 706 K", r"$\Delta^+=4.8$" "\n" r"$s_T/s_{T,th}=0.83$")],
      "H1_pressure", "Same mixing-layer thickness, higher pressure: sharper tip, slower flame")
# H2 squishing at 1 atm
three([("C1", "wide layer (C1)", r"$\Delta^+=32$" "\n" r"$s_T/s_{T,th}=1.07$"),
       ("C4", "squished (C4)", r"$\Delta^+=4.9$" "\n" r"$s_T/s_{T,th}=0.93$"),
       ("C6", "strongly squished (C6)", r"$\Delta^+=1.9$" "\n" r"$s_T/s_{T,th}=0.70$")],
      "H2_squishing", "Squishing the mixing layer at 1 atm: the leading edge sharpens and slows down",
      xw=(-14, 14), yw=(-8, 26))
# H4 inlet
fig, ax = plt.subplots(1, 2, figsize=(8.4, 5.4), sharey=False)
pc = flame_panel(ax[0], "A6P-tanh", "tanh layer, tip 3 mm from inlet", r"$\Delta^+=5.0$, $\varepsilon^*=1.20$" "\n"
                 r"$s_T/s_{T,th}=0.83$", xw=(-12, 12), yw=(-11, 20), show_inlet=True)
flame_panel(ax[1], "A6P-erf", "erf layer, tip far from inlet", r"$\Delta^+=4.9$, $\varepsilon^*=1.62$" "\n"
            r"$s_T/s_{T,th}=0.93$", xw=(-12, 12), yw=(-11, 20))
ax[0].set_ylabel(r"$(y-y_{tip})/\delta_L$")
cb = fig.colorbar(pc, ax=ax, shrink=0.8, pad=0.015); cb.set_label("heat release / max")
fig.canvas.draw()
fig.suptitle("Near the inlet the tip is sharper than its mixing layer predicts, and the flame slower", fontsize=12,
             y=ax[1].get_position().y1 + 0.10)
fig.savefig(os.path.join(OUT, "H4_inlet.png")); plt.close(fig)
print("wrote H4_inlet", flush=True)
# H3 transport (cold runs)
fig = plt.figure(figsize=(11, 5.0))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15], wspace=0.3)
for k, (case, title) in enumerate((("COLDF-C4", "mixture-averaged transport"), ("COLD-C4-SC", "constant Schmidt numbers"))):
    a = fig.add_subplot(gs[k])
    X, Y, F, m = field(case, ["mixture_fraction"], xw=(-12, 12), yw=(0, 100), lev=2, cold=True)
    zf = F["mixture_fraction"]
    a.pcolormesh(X, Y, np.abs(np.gradient(zf, X, axis=1)), cmap="magma", vmin=0, vmax=0.25, shading="auto",
                 rasterized=True)
    lv = [0.01, ZST, 0.1, 0.3, 0.5, 0.7, 0.9]
    a.contour(X, Y, zf, levels=lv, colors="w", linewidths=0.6)
    a.set_xlim(-12, 12); a.set_title(title, fontsize=11); a.set_xlabel(r"$x/\delta_L$")
    if k == 0:
        a.set_ylabel(r"$y/\delta_L$ (from the inlet)")
a = fig.add_subplot(gs[2])
for case, lab, c in (("COLDF-C4", "mixture-averaged", "C0"), ("COLD-C4-SC", "constant Schmidt", "C3")):
    d = json.load(open(os.path.join(PELE, case, "delta_profile.json")))
    y = np.array(d["y"]) / d["deltaL"]; D = np.array(d["Delta"]) / d["deltaL"]
    a.plot(y, D, "o-", color=c, ms=3, label=lab)
a.axvline(18.0e-3 / 317e-6, color="0.5", ls=":", lw=0.8); a.text(18.0e-3 / 317e-6 + 1, 3.6, "tip height (C4)", fontsize=9, color="0.4")
a.set_xlabel(r"$y/\delta_L$"); a.set_ylabel(r"$\Delta^+=\Delta/\delta_L$"); a.legend(frameon=False, fontsize=9)
a.set_title("layer thickness", fontsize=11)
fig.suptitle("The same inlet layer grows 2.4 times faster with constant Schmidt numbers (cold flow, case C4;\n"
             "colour: mixture-fraction gradient, lines: z = 0.01, 0.028 (stoichiometric), 0.1 ... 0.9)", fontsize=12, y=1.04)
fig.savefig(os.path.join(OUT, "H3_transport.png")); plt.close(fig)
print("wrote H3_transport", flush=True)
