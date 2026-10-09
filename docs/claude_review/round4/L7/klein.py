"""NUMERICAL: Klein-symmetric ansatz for rational curves of class beta in Mbar_{0,7}.

Points: anchors 0, 1, inf and movers a, b, c, d on P^1_t, with
  a = A(v), v = t + 1/t   (deg A = 4)        tau_a: t -> -t  moves only a
  b = B(y), y = t^2        (deg B = 4)        tau_b: t -> 1/t moves only b
  c = C(z), d = D(z), z = t^2 + t^-2 = v^2 - 2 = y + 1/y   (deg 2)
Anchor fibre pattern at each anchor value x: four Klein orbits over z_1..z_4 with
  C = x on z_1, z_2;  D = x on z_2, z_3;  A = x at v_i over z_i (i=1..4);  B = x at y_i over z_i.
Non-anchor events: O_cd (c=d=a=b at one point), Q_abc, Q_abd as in the write-up.
"""
import numpy as np, sys

X = [0, 1, 2]  # 0, 1, inf

def monic(roots):
    return np.poly(roots)  # highest first

def unpack(u):
    v = u[0:12].reshape(3, 4); y = u[12:24].reshape(3, 4)
    kA, kB, kC, kD = u[24:28]; vcd, ycd = u[28], u[29]
    return v, y, kA, kB, kC, kD, vcd, ycd

def maps(u):
    v, y, kA, kB, kC, kD, vcd, ycd = unpack(u)
    z = v**2 - 2
    PA = [monic(v[x]) for x in X]; PB = [monic(y[x]) for x in X]
    QC = [monic(z[x, [0, 1]]) for x in X]; QD = [monic(z[x, [1, 2]]) for x in X]
    return v, y, z, PA, PB, QC, QD, kA, kB, kC, kD, vcd, ycd

def ev(p, x):
    return np.polyval(p, x)

def residual(u):
    v, y, z, PA, PB, QC, QD, kA, kB, kC, kD, vcd, ycd = maps(u)
    R = []
    # z-compatibility
    R += list((v**2 - 2 - y - 1/y).ravel())
    # dependencies:  k P0 - Pinf - (k-1) P1 == 0 (drop leading coefficient)
    for (P, k) in [(PA, kA), (PB, kB)]:
        dep = k*P[0] - P[2] - (k-1)*P[1]
        R += list(dep[1:])
    for (Q, k) in [(QC, kC), (QD, kD)]:
        dep = k*Q[0] - Q[2] - (k-1)*Q[1]
        R += list(dep[1:])
    # functions
    A = lambda w: kA*ev(PA[0], w)/ev(PA[2], w)
    B = lambda w: kB*ev(PB[0], w)/ev(PB[2], w)
    C = lambda w: kC*ev(QC[0], w)/ev(QC[2], w)
    D = lambda w: kD*ev(QD[0], w)/ev(QD[2], w)
    # z_cd: remaining root of kC*QC0*QD_inf - kD*QD0*QC_inf (quartic)
    N = kC*np.polymul(QC[0], QD[2]) - kD*np.polymul(QD[0], QC[2])
    sroots = -N[1]/N[0]
    zcd = sroots - z[0, 1] - z[1, 1] - z[2, 1]
    gam = C(zcd)
    R += [vcd**2 - 2 - zcd, ycd + 1/ycd - zcd, A(vcd) - gam, B(ycd) - gam]
    # remaining root of A(v) = C(v^2-2): N_AC = kA*PA0*QCinf(v^2-2) - kC*QC0(v^2-2)*PAinf  (deg 8)
    sq = np.array([1, 0, -2])
    def comp(Q):  # Q(v^2-2) as poly in v
        out = np.zeros(1)
        for co in Q:
            out = np.polyadd(np.polymul(out, sq), [co])
        return out
    def compy(Q):  # Q(y+1/y)*y^deg as poly in y
        deg = len(Q)-1
        out = np.zeros(1)
        # sum co_k (y+1/y)^(deg-k) * y^deg  -> co_k (y^2+1)^(deg-k) y^k
        for k, co in enumerate(Q):
            term = np.array([1.0])
            for _ in range(deg-k):
                term = np.polymul(term, [1, 0, 1])
            term = np.polymul(term, np.concatenate([[1.0], np.zeros(k)]))
            out = np.polyadd(out, co*term)
        return out
    NAC = kA*np.polymul(PA[0], comp(QC[2])) - kC*np.polymul(comp(QC[0]), PA[2])
    NAD = kA*np.polymul(PA[0], comp(QD[2])) - kD*np.polymul(comp(QD[0]), PA[2])
    NBC = kB*np.polymul(PB[0], compy(QC[2])) - kC*np.polymul(compy(QC[0]), PB[2])
    NBD = kB*np.polymul(PB[0], compy(QD[2])) - kD*np.polymul(compy(QD[0]), PB[2])
    def remaining(Npoly, known):
        return -Npoly[1]/Npoly[0] - sum(known)
    # known roots of NAC: v over z_1, z_2 for each anchor, and vcd
    vac = remaining(NAC, [v[x, 0] for x in X] + [v[x, 1] for x in X] + [vcd])
    vad = remaining(NAD, [v[x, 1] for x in X] + [v[x, 2] for x in X] + [vcd])
    ybc = remaining(NBC, [y[x, 0] for x in X] + [y[x, 1] for x in X] + [ycd])
    ybd = remaining(NBD, [y[x, 1] for x in X] + [y[x, 2] for x in X] + [ycd])
    R += [vac**2 - 2 - ybc - 1/ybc, vad**2 - 2 - ybd - 1/ybd]
    return np.array(R), dict(vac=vac, vad=vad, ybc=ybc, ybd=ybd, zcd=zcd, gam=gam,
                             NAC=NAC, NAD=NAD, NBC=NBC, NBD=NBD)

def jac(u, h=1e-7):
    F0, _ = residual(u)
    J = np.zeros((len(F0), len(u)), complex)
    for k in range(len(u)):
        du = np.zeros(len(u), complex); du[k] = h
        J[:, k] = (residual(u+du)[0] - residual(u-du)[0])/(2*h)
    return F0, J

def newton(u, iters=200):
    for it in range(iters):
        F, J = jac(u)
        nF = np.linalg.norm(F)
        if not np.isfinite(nF):
            return None, it
        if nF < 1e-11:
            # polish
            for _ in range(3):
                F, J = jac(u); u = u + np.linalg.solve(J, -F)
            return u, it
        try:
            du = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError:
            return None, it
        st = 1.0
        while st > 1e-8:
            u1 = u + st*du
            F1, _ = residual(u1)
            if np.all(np.isfinite(F1)) and np.linalg.norm(F1) < nF*(1-0.05*st):
                break
            st /= 2
        if st <= 1e-8:
            return None, it
        u = u1
    return None, iters

if __name__ == '__main__':
    nst = int(sys.argv[1]); start = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    sols = []
    for seed in range(start, start+nst):
        rng = np.random.default_rng(seed)
        u = rng.normal(size=30) + 1j*rng.normal(size=30)
        u, it = newton(u)
        if u is None:
            continue
        np.save('klein_sol_%d.npy' % seed, u)
        print(seed, 'converged', it, flush=True)
        sols.append(seed)
    print('converged seeds', sols)
