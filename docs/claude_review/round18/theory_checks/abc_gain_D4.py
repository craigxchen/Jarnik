"""Degree-4 unit pencils, A up to S_7, with a rigorous pruning bound.

Interior factors Y have Y.beta >= 1 (Y.beta is a positive integer: beta is in the interior of the
movable cone, L8.md Thm A(a), refereed).  Ignoring interior contributions to nu_T(B) can only
shrink the boundary support count.  Hence
    gain <= cheap := L.beta - #supp(boundary, pair part of B only) - #(interior Q-irreducible factors).
Pencils with cheap <= THRESH have gain <= THRESH and are not evaluated further.
Partial results are valid as reported (progress per A-orbit).
"""
import abc_gain as G

THRESH = -40
D = 4


def cheap(A, C, s):
    PA, PC = G.mono_poly(A), G.mono_poly(C)
    B = PC - s * PA
    if B == 0:
        return None
    c, facs = B.factor()
    nuB = {T: 0 for T in G.VIS}
    ninter = 0
    for (Fk, mult) in facs:
        ispair = False
        if Fk.total_degree() == 1:
            for (i, j) in G.PAIRS:
                if Fk == G.X[i] - G.X[j] or Fk == G.X[j] - G.X[i]:
                    ispair = True
                    for T in G.VIS:
                        if i in T and j in T:
                            nuB[T] += mult
                    break
        if not ispair:
            ninter += 1
    Lb = 31 * D
    supp = 0
    for T in G.VIS:
        a, cc = G.nu_mono(A, T), G.nu_mono(C, T)
        mu = min(a, cc) if len(T) >= 3 else 0
        if len(T) >= 3:
            Lb -= mu
        if max(a, cc, nuB[T]) > mu:
            supp += 1
    return Lb - supp - ninter


if __name__ == '__main__':
    allm = list(G.monos(D))
    reps = {}
    for m in allm:
        reps.setdefault(G.canon(m), m)
    print(f"D=4: {len(reps)} A-orbits, {len(allm)} monomials; exact gain when cheap bound > {THRESH}", flush=True)
    best = None
    n = nexact = 0
    hist = {}
    cmax = None
    for k, (key, A) in enumerate(reps.items()):
        for C in allm:
            if set(A) & set(C):
                continue
            for s in (1, -1):
                cb = cheap(A, C, s)
                if cb is None:
                    continue
                n += 1
                cmax = cb if cmax is None else max(cmax, cb)
                if cb <= THRESH:
                    continue
                g = G.gain(A, C, s)
                nexact += 1
                hist[g[0]] = hist.get(g[0], 0) + 1
                if best is None or g[0] > best[0]:
                    best = (g[0], g, A, C, s)
        print(f"orbit {k+1}/{len(reps)} A={A}: pencils {n}, exact evaluations {nexact}, "
              f"max cheap bound {cmax}, max exact gain {best[0] if best else None} {best[1:] if best else ''}",
              flush=True)
    print("gains among exactly evaluated:", sorted(hist.items(), reverse=True)[:8], flush=True)
