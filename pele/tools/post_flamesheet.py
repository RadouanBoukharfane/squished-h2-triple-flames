"""PeleLMeX 1D flame speed from fextract lines (lines/plt*.dat: y, temp, y_velocity along the centreline):
y_f = position of T = (Tu + Tb)/2, sL = u_inlet - dy_f/dt fitted on the second half of the saves; also
deltaL = (Tb - Tu)/max dT/dy and the drift between the two quarters of the fit window.
Usage: python post_flamesheet.py LINES_DIR [SL_CANTERA]"""
import glob, sys
import numpy as np

t, yf, uin, dl = [], [], [], []
for fn in sorted(glob.glob(sys.argv[1] + "/plt*.dat")):
    tt = float(open(fn).readlines()[1].split("=")[1])
    y, T, v = np.loadtxt(fn, unpack=True)
    Tm = 0.5 * (T[0] + T[-1])
    k = np.argmax(T > Tm)
    t.append(tt); uin.append(v[0]); dl.append((T[-1] - T[0]) / np.max(np.gradient(T, y)))
    yf.append(y[k - 1] + (Tm - T[k - 1]) * (y[k] - y[k - 1]) / (T[k] - T[k - 1]))
t, yf, uin, dl = map(np.array, (t, yf, uin, dl))
for a, b in zip(t, yf):
    print("  t = %.3e s  y_f = %.4f mm" % (a, 1e3 * b))
w = slice(len(t) // 2, None)
sL = np.mean(uin[w]) - np.polyfit(t[w], yf[w], 1)[0]
q = len(t) - len(t) // 2
a, b = slice(len(t) // 2, len(t) // 2 + q // 2), slice(len(t) // 2 + q // 2, None)
dr = [np.mean(uin[s]) - np.polyfit(t[s], yf[s], 1)[0] for s in (a, b)]
msg = "sL_Pele = %.4f m/s  deltaL = %.1f um  drift = %+.2f %%" % (sL, 1e6 * np.mean(dl[w]), 100 * (dr[1] - dr[0]) / sL)
if len(sys.argv) > 2:
    msg += "  vs Cantera %.4f: %+.2f %%" % (float(sys.argv[2]), 100 * (sL / float(sys.argv[2]) - 1))
print(msg)
