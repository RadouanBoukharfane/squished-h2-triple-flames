"""Cold companions of the runs auto_cold.py does not handle (29/10/2026): B5-MA-D20, B2-MA-D20 and MESH-p10D5b-L3.
Same construction as auto_cold.py (make_pele_case.py --cold, same Delta0+ and height, V = final inlet velocity,
mixture-averaged transport + regrid_interp_method 0 / floor_species 1, 4 nodes if n_cell 960). For a MESH case, every
amr.*max_level of the cold run is capped at MAXL, so that its layer resolution follows the reactive case (the cold
base grid is twice the reactive one: cold level l = reactive level l + 1). Every 10 min; once the Slurm job has left
the queue: stop time reached -> build and submit COLDF-<case>; otherwise log FAILED. Log: PELE/auto_cold.log.
Usage: setsid nohup python tools/auto_cold_manual.py >/dev/null 2>&1 < /dev/null &"""
import os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))          # PELE/
# case: (Delta0+, y_target mm (1-atm eq.), p bar, Tu K, cap on max_level or None)
CASES = {"B5-MA-D20": (19, 18.7, 5, 475, None),
         "B2-MA-D20": (19, 18.7, 10, 300, None),
         "MESH-p10D5b-L3": (2, 13.3, 10, 579, 2)}                          # reactive max_level 3 (base 4)
LOG = os.path.join(HERE, "auto_cold.log")


def log(msg):
    open(LOG, "a").write(time.strftime("%d/%m %H:%M ") + "[manual] " + msg + "\n")


def queued(case):
    out = subprocess.run(["squeue", "-h", "-u", os.environ["USER"], "-o", "%j"], capture_output=True, text=True).stdout
    return ("pele_" + case) in out.split()


done = set()
log("start: " + " ".join(CASES))
while len(done) < len(CASES):
    for case, (d0, yt, p, tu, maxl) in CASES.items():
        if case in done or queued(case):
            continue
        done.add(case)
        cold = "COLDF-" + case
        if os.path.exists(os.path.join(HERE, cold)):
            log("%s exists: nothing done" % cold); continue
        inp = open(os.path.join(HERE, case, "input.case")).read()
        stop = float(re.search(r"^amr.stop_time\s*=\s*(\S+)", inp, re.M).group(1))
        rows = [l.split() for l in open(os.path.join(HERE, case, "AC_History.dat")) if l.strip()]
        t, V = float(rows[-1][1]), float(rows[-1][2])
        if t < 0.999 * stop:
            log("%s FAILED or stopped early (t = %.4e of %.4e): no cold run" % (case, t, stop)); continue
        cmd = [sys.executable, os.path.join(HERE, "tools", "make_pele_case.py"), cold, repr(d0), repr(yt), "%.5f" % V,
               "--cold", "--hours", "6", "--p", repr(p), "--tu", repr(tu)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            log("%s: make_pele_case failed: %s" % (cold, r.stderr[-300:])); continue
        fi = os.path.join(HERE, cold, "input.case")
        s = open(fi).read()
        if maxl is not None:
            s = re.sub(r"^(amr\.(?:\w+\.)?max_level\s*=\s*)(\d+)", lambda m: m.group(1) + str(min(int(m.group(2)), maxl)),
                       s, flags=re.M)
            s = re.sub(r"^(amr\.(?:ref_ratio|n_error_buf)\s*=\s*)(.*)$", lambda m: m.group(1) + " ".join(["2"] * maxl),
                       s, flags=re.M)
            s += "# MESH cold run: every max_level capped at %d (layer resolution of %s)\n" % (maxl, case)
        s += "peleLM.regrid_interp_method = 0\npeleLM.floor_species = 1\n"
        open(fi, "w").write(s)
        sb = os.path.join(HERE, cold, "run.sbatch")
        if re.search(r"^amr.n_cell\s*=\s*960", s, re.M):
            b = open(sb).read(); open(sb, "w").write(b.replace("--nodes=2", "--nodes=4"))
        j = subprocess.run(["sbatch", "--parsable", "run.sbatch"], cwd=os.path.join(HERE, cold), capture_output=True,
                           text=True)
        log("%s: V = %.4f m/s (t = %.4e) -> %s submitted, job %s" % (case, V, t, cold, j.stdout.strip() or j.stderr.strip()))
    time.sleep(600)
log("all cases handled")
