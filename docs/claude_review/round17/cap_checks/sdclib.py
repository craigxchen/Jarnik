# shared sdc_q construction (independent re-implementation of the Paley-type skew double core, as in
# round16/referee_barrier_checks/sdc_sym.py). coordinates x in F_q (base row inf omitted).
def chi(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1
def patterns(q):
    """P_a = {a} u (a+NQR) ; Pp_a = {a} u (a+QR)  (q = 3 mod 4)."""
    P = [frozenset([a] + [x for x in range(q) if chi(a - x, q) == 1]) for a in range(q)]
    Pp = [frozenset([a] + [x for x in range(q) if chi(a - x, q) == -1]) for a in range(q)]
    return P, Pp
