# Recompute 60 random classes sigma_*gamma_A of lamA_bound.py with a second, independent V_A-member
# (evidence that the counts are generic, i.e. that the test curves avoid Sing(F_A)).
import ast, random
import lamA_bound as B, lamA_exact as E
d = {}
for k in range(3):
    d.update(dict(ast.literal_eval(open('lamA_bound_part%d_of3.txt' % k).read())['vals']))
rng = random.Random(9)
for s in range(5001, 5100):
    try:
        ok, mob, x = B.member_mob(s)
    except Exception: continue
    if ok: break
print('second member through', [str(c) for c in x], 'class verified', ok)
same = diff = 0
for sig in rng.sample(sorted(d), 60):
    v = E.lam_from_br(*E.leaf_brackets([mob[sig.index(a)] for a in range(7)]))
    if v == d[sig]: same += 1
    else: diff += 1; print('differs', sig, d[sig], v)
print('same', same, 'differ', diff)
