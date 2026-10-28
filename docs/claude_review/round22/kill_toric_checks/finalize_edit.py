p = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round22/kill_toric.md'
s = open(p).read()
old = "  * *`k = 5`:* see §8 for the result. It does not affect any verdict, because nothing depends on it."
new = ("  * *`k = 5`* (`pte_own2.py 5 80 11 3000` → `pte_k5b.out`): 44/80 converged and **0 were nondegenerate**. In every\n"
       "    converged solution some roots coincide: the largest minimal root distance is `2·10^{-12}`. This **matches** the\n"
       "    author's finding in kind. A longer run (`pte_own.py 5 400 7`, `max_nfev = 20000`) was stopped after about 30 CPU\n"
       "    minutes without output and is not counted. Evidence only. Nothing depends on it.")
assert old in s
s = s.replace(old, new)
o2 = "| `pte_own.py` | `pte_k5.out` (k=5); k=4 calibration printed in the session (100/100 converged, 24 nondegenerate) |"
assert o2 in s
s = s.replace(o2, "| `pte_own.py`, `pte_own2.py` (same, flushed, `max_nfev` argument) | `pte_k5b.out` (k=5: 44/80 converged, 0 nondegenerate); k=4 calibration printed in the session (100/100 converged, 24 nondegenerate) |")
s = s.replace("Status: IN PROGRESS\n", "Status: complete (killer report, round 22). Nothing here is refereed.\n")
o3 = "7. **(evidence)** No nondegenerate real solution of the minimal-degree `k = 5` PTE system has been found (§3, and the\n   result recorded below)."
assert o3 in s
s = s.replace(o3, "7. **(evidence)** No nondegenerate real solution of the minimal-degree `k = 5` PTE system has been found, in the\n   author's search or in mine (§3); the `k = 4` calibration finds solutions readily.")
open(p, 'w').write(s)
print("ok", s.count("(pending)"))
