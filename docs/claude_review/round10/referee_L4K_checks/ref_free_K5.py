from ref_free import run
ps = [5, 13, 17, 29]
f = run(ps, 5, 420)
print(f"K=5 r=4 first ps={ps} q<=420: q*={max(f)} #forced={len(f)} forced={f}", flush=True)
