"""Referee: independent recomputation of Cor 2.3 (n0, log mu bound), Cor 5.2 thresholds M_0, and tau of §8.1."""
import mpmath as mp, json, math
mp.mp.dps = 40
out = {}
# --- tau for the §8.1 prime set
ps = [13,17,29,37,41,53,61,73]; split=[5,13,17,29,37,41,53,61,73]
qs = [split[split.index(p)-1] for p in ps]
ratios = {f'{p}/{q}': p/q for p,q in zip(ps,qs)}
tau2 = min(ratios.values()); out['sec8.1'] = {'ratios': ratios, 'tau^2=min p/q': tau2, 'tau': math.sqrt(tau2), 'claimed tau^2': 73/61}
# --- Cor 2.3: n0 = least integer with (c2 r/2)(tau/(tau-1)) tau^{-n0/2} f(z)(1+2rz/(1-z^2)) <= 1/2
c2 = mp.pi*(1+mp.sqrt(2))
def cor23(r, P, tau):
    tau = mp.mpf(tau); z = tau**-0.5
    A = (c2*r/2)*(tau/(tau-1))*((1+z)/(1-z))**r*(1+2*r*z/(1-z*z))
    n0 = 1
    # tau^{-n0/2} A <= 1/2  <=> n0 >= 2 log(2A)/log tau
    n0 = int(mp.ceil(2*mp.log(2*A)/mp.log(tau)))
    while (tau**(-mp.mpf(n0-1)/2))*A <= 0.5: n0 -= 1
    while (tau**(-mp.mpf(n0)/2))*A > 0.5: n0 += 1
    eps = min(mp.pi/5, (tau-1)*mp.mpf(P)**(-mp.mpf(n0-1)/2)/(mp.sqrt(2)*n0))
    logmu = mp.log(0.5) + r*mp.log(4*eps/mp.pi)
    asym = mp.log((1+z)/(1-z))/mp.log(tau)
    return n0, float(logmu), float(asym)
rows = []
for (r,P,tau) in [(8,73,math.sqrt(73/61)),(8,73,math.sqrt(41/37)),(8,1000,2),(100,10**4,2),(2500,10**6,2)]:
    n0, lm, a = cor23(r,P,tau)
    rows.append({'r':r,'P':P,'tau':tau,'n0':n0,'log_mu_lower':lm,'asym_const':a})
out['cor23'] = rows
# Poisson comparison for r=2500 primes ~1e6: -1.1(Pi-1), Pi=(1+2/(sqrt P-1))^r
Pi = (1+2/(mp.sqrt(10**6)-1))**2500; out['poisson_r2500'] = float(-1.1*(Pi-1))
# --- Cor 5.2: last zero of f(M) = c M^{d/2} - log(1+2M^2/kappa)
kappa = (4/mp.pi)*(1+1/mp.sqrt(2))
def M0(c,d):
    f = lambda L: c*mp.e**(L*d/2) - mp.log(1+2*mp.e**(2*L)/kappa)   # L = log M
    # scan from large L downward to find last sign change
    Ls = [mp.mpf(k)/10 for k in range(0, 1500)]
    vals = [f(L) for L in Ls]
    last = None
    for i in range(len(Ls)-1):
        if vals[i] <= 0 < vals[i+1]: last = i
    assert all(v>0 for v in vals[last+1:])
    L = mp.findroot(f, (Ls[last], Ls[last+1]), solver='bisect')
    return float(mp.e**L)
out['M0'] = {f'c={c},delta={d}': M0(c,d) for (c,d) in [(0.25,1.0),(1.0,0.5),(0.25,0.5),(0.5,0.25),(0.25,0.25)]}
print(json.dumps(out, indent=1))
json.dump(out, open('constants_out.json','w'), indent=1)
