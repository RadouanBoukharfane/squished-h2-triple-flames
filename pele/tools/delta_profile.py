"""Delta(y) = 1 / max dz/dx of a PeleLMeX cold run (last plotfile), z = mixture_fraction, via fextract
slices in x at several heights; prints Delta+ = Delta/deltaL and fits Delta+^2 = a + b Y+.
Usage (foss/2024a loaded for fextract): python delta_profile.py CASE_DIR [N_HEIGHTS=12]"""
import glob, json, os, re, subprocess, sys
import numpy as np

FX = os.path.join(os.environ["PELE_HOME"], "Submodules/PelePhysics/Submodules/amrex/Tools/Plotfile/fextract.gnu.ex")
case = sys.argv[1]
n = int(sys.argv[2]) if len(sys.argv) > 2 else 12
inp = open(os.path.join(case, "input.case")).read()
get = lambda k: float(re.search(r"^%s\s*=\s*(\S+)" % re.escape(k), inp, re.M).group(1))
ly = float(re.search(r"^geometry.prob_hi\s*=\s*\S+\s+(\S+)", inp, re.M).group(1))
dl0 = get("prob.delta0") / float(re.search(r"Delta0\+ = ([0-9.]+)", inp).group(1))     # deltaL from the header
plt = sorted(p for p in glob.glob(os.path.join(case, "plt?????")))[-1]
rows = []
for y in np.linspace(0.02, 0.98, n) * ly:
    out = os.path.join(case, "slice_y%.6e.dat" % y)
    subprocess.run([FX, "-d", "0", "-y", "%.8e" % y, "-v", "mixture_fraction", "-s", out, plt], check=True,
                   stdout=subprocess.DEVNULL)
    x, z = np.loadtxt(out, unpack=True)
    rows.append((y, 1.0 / np.max(np.gradient(z, x))))
rows = np.array(rows)
b, a = np.polyfit(rows[:, 0] / dl0, (rows[:, 1] / dl0) ** 2, 1)
print("%s (%s): deltaL = %.2f um" % (case, os.path.basename(plt), dl0 * 1e6))
for y, d in rows:
    print("  y = %8.4f mm  Y+ = %6.1f  Delta+ = %6.2f" % (1e3 * y, y / dl0, d / dl0))
print("  fit: Delta+^2 = %.2f + %.3f Y+   (Delta0+ fit = %.2f)" % (a, b, np.sqrt(max(a, 0))))
json.dump(dict(case=case, plt=plt, deltaL=dl0, y=list(rows[:, 0]), Delta=list(rows[:, 1]), a=a, b=b),
          open(os.path.join(case, "delta_profile.json"), "w"), indent=1)
