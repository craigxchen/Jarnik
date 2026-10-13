from ref_free import run
from rlib import split_primes_upto
sp = split_primes_upto(200)
for which, ps in (("first", sp[:4]), ("last", sp[-4:])):
    f = run(ps, 4, 700)
    print(f"K=4 r=4 {which} ps={ps} q<=700: q*={max(f)} #forced={len(f)} forced={f}", flush=True)
