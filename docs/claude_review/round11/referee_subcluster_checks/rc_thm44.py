# Independent exact test of Theorem 4.4: reg(S) <= Omega'(S) = E_eps h_S(<eps,.>)
import sys, random
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round11/referee_subcluster_checks')
from rc_lib import *
random.seed(20261013)
worst = {}
viol = 0
count = 0
# random sets
for dim, box, sizes, trials in [(1,6,(2,7),100),(2,4,(2,20),400),(3,3,(2,25),300),(4,2,(2,20),150)]:
    pts = list(product(range(box+1), repeat=dim))
    for t in range(trials):
        s = random.randint(*sizes)
        s = min(s, len(pts))
        S = random.sample(pts, s)
        rg = reg(S); om = omega_signed(S); mw = MW(S)
        count += 1
        if rg > om: viol += 1; print("VIOLATION", dim, S, rg, om)
        if rg > mw: print("THM B VIOLATION", dim, S)
        ratio = Fraction(rg) / om if om > 0 else 0
        if ratio > worst.get(dim, (0,None))[0]:
            worst[dim] = (ratio, S, rg, om)
print("random sets tested:", count, "violations:", viol)
for d in sorted(worst):
    r, S, rg, om = worst[d]
    print("dim", d, "max reg/Omega' =", float(r), "example size", len(S), "reg", rg, "Omega'", om)
# adversarial hill climbing in dim 2 and 3: maximise reg - Omega'
def score(S):
    return reg(S) - omega_signed(S)
best_overall = None
for dim, box in [(2,5),(3,3)]:
    pts = list(product(range(box+1), repeat=dim))
    for restart in range(25):
        S = set(random.sample(pts, random.randint(3, 10)))
        sc = score(list(S))
        for it in range(150):
            T = set(S)
            mv = random.random()
            if mv < 0.4 and len(T) > 2:
                T.remove(random.choice(list(T)))
            elif mv < 0.8:
                T.add(random.choice(pts))
            else:
                if len(T) > 2:
                    T.remove(random.choice(list(T))); T.add(random.choice(pts))
            if len(T) < 2: continue
            s2 = score(list(T))
            if s2 >= sc:
                S, sc = T, s2
        if best_overall is None or sc > best_overall[0]:
            best_overall = (sc, dim, sorted(S))
        if sc > 0: print("VIOLATION (hill)", dim, sorted(S), sc)
print("hill-climb best reg - Omega':", float(best_overall[0]), "dim", best_overall[1], "set", best_overall[2])
# equality cases claimed
for dim, R in [(2,2),(2,3),(3,2),(3,3)]:
    S = [y for y in product(range(-R,R+1), repeat=dim) if sum(map(abs,y)) <= R]
    print("l1-ball dim", dim, "R", R, "reg", reg(S), "Omega'", omega_signed(S), "MW", MW(S))
for a in [(1,2),(2,3),(1,1,2)]:
    S = list(product(*[range(x+1) for x in a]))
    print("box", a, "reg", reg(S), "Omega'", omega_signed(S), "sum", sum(a))
