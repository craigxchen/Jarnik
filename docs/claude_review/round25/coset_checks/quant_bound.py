"""Corollary 2.3: explicit lower bound log mu(G(p)) >= log(1/2) + r log(4 eps/pi) from Lemma 2.1 (full torus),
eps = (tau-1) P^{-(n0-1)/2} / (sqrt2 * n0), n0 minimal with
 F_bound(n0) = (c2 r/2) (tau/(tau-1)) tau^{-n0/2} ((1+z1)/(1-z1))^r (1 + 2 r z1/(1-z1^2)) <= 1/2, z1 = tau^{-1/2}.
Compares with Theorem A (exp(-kappa Lambda)) and the Poisson extrapolation exp(-1.1 (Pi-1)) of fakes.md 5.1.
All logs natural."""
import math, json
c2 = math.pi * (1 + math.sqrt(2))
def bound(r, P, tau):
    z1 = tau ** -0.5
    logF0 = math.log(c2 * r / 2) + math.log(tau / (tau - 1)) + r * math.log((1 + z1) / (1 - z1)) + math.log(1 + 2 * r * z1 / (1 - z1 * z1))
    n0 = max(1, math.ceil(2 * (logF0 - math.log(0.5)) / math.log(tau)))
    logeps = math.log(tau - 1) - (n0 - 1) / 2 * math.log(P) - math.log(math.sqrt(2) * n0)
    return n0, math.log(0.5) + r * (math.log(4 / math.pi) + logeps)
rows = []
for (r, P, tau) in [(8, 73, (73/61) ** 0.5), (8, 1000, 2.0), (100, 10**4, 2.0), (2500, 10**6, 2.0),
                    (10**4, 10**8, 2.0), (10**6, 10**12, 2.0)]:
    n0, lb = bound(r, P, tau)
    Lam = r / math.sqrt(P)                     # primes ~P: Lambda ~ r P^{-1/2} (upper end; primes in [P/2,P] give <= sqrt2 times this)
    logPi = r * math.log(1 + 2 / (math.sqrt(P) - 1))
    rows.append({'r': r, 'P': P, 'tau': round(tau, 4), 'n0': n0, 'log_mu_lower': lb,
                 'thmA_log_upper': -2.1735 * Lam, 'poisson_log': -1.1 * (math.exp(logPi) - 1) if logPi < 700 else -float('inf'),
                 'logPi': logPi})
for x in rows: print(x)
json.dump(rows, open('quant_bound_out.json', 'w'), indent=1)
