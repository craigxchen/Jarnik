# finalize.py -- regenerate rigcompute.md from rigcompute.tmpl.md with the current tables.
import subprocess, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
run = lambda *a: subprocess.run(["python3", *a], capture_output=True, text=True).stdout
t = open("../rigcompute.tmpl.md").read()
t = t.replace("WITNESS_EXAMPLES", run("wit_examples.py"))
t = t.replace("WITNESS_TABLE", run("witness_table.py"))
t = t.replace("JOINT_TABLE", run("make_joint_table.py"))
t = t.replace("THR_TABLE", run("thr_table.py"))
t = t.replace("ONEMOD_TABLE", run("summarize.py"))
open("../rigcompute.md", "w").write(t)
print("ok", t.count("PENDING"))
