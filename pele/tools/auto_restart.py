"""Automatic restart of runs stopped by amr.max_wall_time before amr.stop_time (29/10/2026). Every 10 min, for each
case given: once its Slurm job has left the queue, if the last AC_History.dat time is below stop_time and pele.log
shows no abort, keep pele.log and AC_History.dat (suffix .partN), add amr.restart = <last checkpoint> to input.case
and resubmit run.sbatch; at most 2 restarts per case. Cold runs are not built here. Log: PELE/auto_restart.log.
Usage: setsid nohup python tools/auto_restart.py CASE [CASE ...] &"""
import glob, os, re, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))          # PELE/
LOG = os.path.join(HERE, "auto_restart.log")


def log(msg):
    open(LOG, "a").write(time.strftime("%d/%m %H:%M ") + msg + "\n")


def queued(case):
    out = subprocess.run(["squeue", "-h", "-u", os.environ["USER"], "-o", "%j"], capture_output=True, text=True).stdout
    return ("pele_" + case) in out.split()


todo = {c: 0 for c in sys.argv[1:]}
log("start: " + " ".join(todo))
while todo:
    for case in list(todo):
        if queued(case):
            continue
        d = os.path.join(HERE, case)
        inp = open(os.path.join(d, "input.case")).read()
        stop = float(re.search(r"^amr.stop_time\s*=\s*(\S+)", inp, re.M).group(1))
        rows = [l.split() for l in open(os.path.join(d, "AC_History.dat")) if l.strip()]
        t, V = float(rows[-1][1]), float(rows[-1][2])
        if t >= 0.999 * stop:
            log("%s finished (t = %.4e, V = %.4f)" % (case, t, V)); del todo[case]; continue
        if re.search(r"abort|nan", open(os.path.join(d, "pele.log")).read(), re.I) or todo[case] >= 2:
            log("%s stopped at t = %.4e without wall-time stop or after 2 restarts: left as is" % (case, t))
            del todo[case]; continue
        chk = sorted(glob.glob(os.path.join(d, "chk?????")))[-1]
        todo[case] += 1
        for f in ("pele.log", "AC_History.dat"):
            shutil.copy(os.path.join(d, f), os.path.join(d, "%s.part%d" % (f, todo[case])))
        inp = re.sub(r"^amr.restart\s*=.*\n", "", inp, flags=re.M)
        open(os.path.join(d, "input.case"), "w").write(inp + "amr.restart = %s\n" % os.path.basename(chk))
        j = subprocess.run(["sbatch", "--parsable", "run.sbatch"], cwd=d, capture_output=True, text=True)
        log("%s: t = %.4e of %.4e, V = %.4f -> restart %d from %s, job %s" %
            (case, t, stop, V, todo[case], os.path.basename(chk), j.stdout.strip() or j.stderr.strip()))
    time.sleep(600)
log("all cases handled")
