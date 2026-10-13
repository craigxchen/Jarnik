"""Thm 3.2 / Prop 3.1 on the record clusters (Cmax 8, 5.3, 5.9) and all runs of those N at Cmax 8."""
from tangent import analyse, all_runs

out = []
for N, Cm in [(1176852625, 8), (235370525, 5.3), (537606725, 5.9)]:
    runs = all_runs(N, Cm)
    m = max(len(r) for r in runs)
    for r in runs:
        v = analyse(r, "rec", out)
        if v: print("VIOLATION", N, v)
    print(N, "runs", len(runs), "max M", m)
for r in sorted(out, key=lambda r: -r[2])[:6]:
    print(f"N={r[1]} M={r[2]} C={r[3]:.3f} W={r[4]:.1f} h*={r[5]:.2f} v*={r[6]} C^2h*+2={r[3]**2*r[5]+2:.1f} #exc={r[7]} #big={r[8]} slack_b={r[9]:.2f} viol={r[10]}")
