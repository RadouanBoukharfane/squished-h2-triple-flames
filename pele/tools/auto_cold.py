"""Automatic cold companions of the PeleLMeX constant-Schmidt campaign (28/10/2026). Every 10 min, for each reactive
case of build_sc_campaign.py (plus C4-SC, B1-p10-D5-SC and the mixture-averaged A6P-tanh10): once its Slurm job has
left the queue, check that it reached its stop time; if so, build COLDF-<case> (make_pele_case.py --cold, same inlet
Delta0+ and height, V = final inlet velocity of the active control, same transport and executable) and submit it;
if not, log it as FAILED and build nothing. Log: PELE/auto_cold.log. Stops when every case is handled.
Usage: nohup python tools/auto_cold.py &"""
import os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))          # PELE/
sys.path.insert(0, HERE)
import build_sc_campaign as B                                               # noqa: E402  (CASES, transport)

CASES = {c[0]: c[1:] for c in B.CASES}
CASES["C4-SC"] = (3.2, 18.0, 7.90, 1, 300)
CASES["B1-p10-D5-SC"] = (3, 7.2, 13.79, 10, 579)
CASES["A6P-tanh10"] = (4.11, 10.0, 7.81, 1, 300)                           # mixture-averaged
LOG = os.path.join(HERE, "auto_cold.log")


def log(msg):
    open(LOG, "a").write(time.strftime("%d/%m %H:%M ") + msg + "\n")


def queued(case):
    out = subprocess.run(["squeue", "-h", "-u", os.environ["USER"], "-o", "%j"], capture_output=True, text=True).stdout
    return ("pele_" + case) in out.split()


def last_ac(case):
    rows = [l.split() for l in open(os.path.join(HERE, case, "AC_History.dat")) if l.strip()]
    return float(rows[-1][1]), float(rows[-1][2])


done = set()
log("start, %d cases" % len(CASES))
while len(done) < len(CASES):
    for case, (d0, yt, v0, p, tu) in CASES.items():
        if case in done or queued(case):
            continue
        cold = "COLDF-" + case
        if os.path.exists(os.path.join(HERE, cold)):
            done.add(case); continue
        inp = open(os.path.join(HERE, case, "input.case")).read()
        stop = float(re.search(r"^amr.stop_time\s*=\s*(\S+)", inp, re.M).group(1))
        t, V = last_ac(case)
        if t < 0.999 * stop:
            log("%s FAILED or stopped early (t = %.4e of %.4e): no cold run" % (case, t, stop))
            done.add(case); continue
        cmd = [sys.executable, os.path.join(HERE, "tools", "make_pele_case.py"), cold, repr(d0), repr(yt), "%.5f" % V,
               "--cold", "--hours", "6"]
        if (p, tu) != (1, 300):
            cmd += ["--p", repr(p), "--tu", repr(tu)]
        sc = case.endswith("-SC")
        if sc:
            cmd += ["--exe-dir", "TripleFlame-SC-build"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            log("%s: make_pele_case failed: %s" % (cold, r.stderr[-300:])); done.add(case); continue
        with open(os.path.join(HERE, cold, "input.case"), "a") as f:
            f.write(B.transport(p, tu) if sc else "peleLM.regrid_interp_method = 0\npeleLM.floor_species = 1\n")
        sb = os.path.join(HERE, cold, "run.sbatch")
        if re.search(r"^amr.n_cell\s*=\s*960", open(os.path.join(HERE, cold, "input.case")).read(), re.M):
            s = open(sb).read(); open(sb, "w").write(s.replace("--nodes=2", "--nodes=4"))
        j = subprocess.run(["sbatch", "--parsable", "run.sbatch"], cwd=os.path.join(HERE, cold), capture_output=True, text=True)
        log("%s: V = %.4f m/s (t = %.4e) -> %s submitted, job %s" % (case, V, t, cold, j.stdout.strip() or j.stderr.strip()))
        done.add(case)
    time.sleep(600)
log("all cases handled")
