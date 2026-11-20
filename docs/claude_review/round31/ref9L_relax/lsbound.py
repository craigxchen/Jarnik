import math
# largest b allowed by the pair-valuation (Vandermonde-refined larger sieve) inequality at L = sqrt q,
# for q all of whose prime factors exceed b:  C(b,2)*W/2 >= W*floor((b-1)^2/4) + sum_{j<b} log j!
def bmax(W):
    b=2
    while True:
        nb=b+1
        lhs=math.comb(nb,2)*W/2; rhs=W*((nb-1)**2//4)+sum(math.lgamma(j+1) for j in range(1,nb))
        if lhs<rhs: return b
        b=nb
for W in [778.2,837.8,897.1,183.9]:
    print(W, bmax(W), round(W/(2*math.log(W)),1))
for e in [3,4,6,8,10]:
    W=10.0**e; print("W=1e%d"%e, bmax(W), "ratio to W/logW:", round(bmax(W)/(W/math.log(W)),3))
