"""Refined cold runs for the measured eps* (review of 29/10/2026, point 5): COLDR-<case> = copy of COLDF-<case> (or
COLD-<case>) with the mixing-layer indicator active up to the finest level and over the whole height
(amr.layer.max_level = amr.max_level, Y(H2) jump per cell lowered from 0.05 to THR), everything else unchanged.
In COLDF-* the layer is refined only to level 2 (about 0.25 deltaL), and only near the inlet where the jump is large.
4 nodes when the original cold run ended above 1M cells. Jobs submitted with --nice=500 (behind the campaign).
Usage: python make_refined_cold.py [--submit] [--thr 0.001] [--hours 24] COLDF-C4 COLDF-p20-D3 ..."""
import argparse, os, re, shutil, subprocess

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("cold", nargs="+"); ap.add_argument("--submit", action="store_true")
ap.add_argument("--thr", type=float, default=0.001); ap.add_argument("--hours", type=int, default=24)
a = ap.parse_args()
for src in a.cold:
    case = re.sub(r"^COLDF?-", "", src)
    dst = "COLDR-" + case
    s, d = os.path.join(HERE, src), os.path.join(HERE, dst)
    if os.path.exists(d):
        print(dst, "exists: skipped"); continue
    os.makedirs(d)
    for f in os.listdir(s):
        if f.endswith(".dat") and not f.startswith(("slice_", "AC_")):
            shutil.copy(os.path.join(s, f), d)
    exe = [f for f in os.listdir(s) if f.endswith(".ex")][0]
    os.symlink(os.path.realpath(os.path.join(s, exe)), os.path.join(d, exe))
    inp = open(os.path.join(s, "input.case")).read()
    inp = re.sub(r"^amr.restart\s*=.*\n", "", inp, flags=re.M)
    # 30/10/2026: max_wall_time of the cold run (3.5 or 5.5 h) aligned on the Slurm limit, it had stopped 6 COLDR early
    inp = re.sub(r"^amr.max_wall_time\s*=.*$", "amr.max_wall_time = %.1f" % (a.hours - 0.5), inp, flags=re.M)
    ml = re.search(r"^amr.max_level\s*=\s*(\d+)", inp, re.M).group(1)
    inp = re.sub(r"^amr.layer.max_level\s*=.*$", "amr.layer.max_level = " + ml, inp, flags=re.M)
    inp = re.sub(r"^amr.layer.adjacent_difference_greater\s*=.*$",
                 "amr.layer.adjacent_difference_greater = %g" % a.thr, inp, flags=re.M)
    inp = ("# Refined cold run %s (PELE/tools/make_refined_cold.py, 29/10/2026): copy of %s with the layer indicator up to\n"
           "# level %s over the whole height (jump %g instead of 0.05), for eps* measured on the cold profile.\n"
           % (dst, src, ml, a.thr)) + inp
    open(os.path.join(d, "input.case"), "w").write(inp)
    sb = open(os.path.join(s, "run.sbatch")).read().replace("pele_" + src, "pele_" + dst)
    sb = re.sub(r"--time=\S+", "--time=%02d:00:00" % a.hours, sb)
    nc = re.findall(r"Ncells (\d+)", open(os.path.join(s, "pele.log")).read())
    if nc and int(nc[-1]) > 1000000:
        sb = re.sub(r"--nodes=\d+", "--nodes=4", sb)
    open(os.path.join(d, "run.sbatch"), "w").write(sb)
    msg = dst + " built (max_level %s)" % ml
    if a.submit:
        j = subprocess.run(["sbatch", "--parsable", "--nice=500", "run.sbatch"], cwd=d, capture_output=True, text=True)
        msg += ", job " + (j.stdout.strip() or j.stderr.strip())
    print(msg)
