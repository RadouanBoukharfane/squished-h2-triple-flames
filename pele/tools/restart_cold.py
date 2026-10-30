"""Restart of cold runs (COLDF-*, COLDR-*) stopped by amr.max_wall_time before amr.stop_time (30/10/2026).
The cold runs inherit amr.max_wall_time from make_pele_case.py (3.5 or 5.5 h) even when run.sbatch asks for 24 h, and
no watcher restarts them. For each run: writes run_restartN.sbatch (copy of run.sbatch, --time HOURS, job name -rN),
which at start keeps pele.log as pele.log.partN, sets amr.restart = last checkpoint and amr.max_wall_time = HOURS - 0.5
in input.case, then runs. --after JOBID: submitted with afterany on a job still running. Jobs submitted with --nice=500.
Usage: python restart_cold.py [--submit] [--hours 24] [--after JOBID] COLDR-S6-erf ..."""
import argparse, glob, os, re, subprocess

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+"); ap.add_argument("--submit", action="store_true")
ap.add_argument("--hours", type=int, default=24); ap.add_argument("--after")
a = ap.parse_args()
for run in a.runs:
    d = os.path.join(HERE, run)
    n = len(glob.glob(os.path.join(d, "run_restart*.sbatch"))) + 1
    sb = open(os.path.join(d, "run.sbatch")).read()
    sb = re.sub(r"--job-name=\S+", "--job-name=pele_%s-r%d" % (run, n), sb)
    sb = re.sub(r"--time=\S+", "--time=%02d:00:00" % a.hours, sb)
    head, tail = sb.split("mpirun", 1)
    head += ("# restart %d (tools/restart_cold.py): last checkpoint, max_wall_time aligned on the Slurm time limit\n"
             "chk=$(ls -d chk????? | sort | tail -1); echo \"restart from $chk\"\n"
             "cp pele.log pele.log.part%d\n"
             "sed -i -e '/^amr.restart *=/d' -e 's/^amr.max_wall_time *=.*/amr.max_wall_time = %.1f/' input.case\n"
             "echo \"amr.restart = $chk\" >> input.case\n" % (n, n, a.hours - 0.5))
    f = "run_restart%d.sbatch" % n
    open(os.path.join(d, f), "w").write(head + "mpirun" + tail)
    msg = "%s: %s written" % (run, f)
    if a.submit:
        cmd = ["sbatch", "--parsable", "--nice=500"] + (["--dependency=afterany:" + a.after] if a.after else []) + [f]
        j = subprocess.run(cmd, cwd=d, capture_output=True, text=True)
        msg += ", job " + (j.stdout.strip() or j.stderr.strip())
    print(msg)
