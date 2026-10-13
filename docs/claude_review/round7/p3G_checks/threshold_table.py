import json, collections, math
rows=[json.loads(l) for l in open('threshold_out.jsonl')]
agg=collections.defaultdict(list)
for d in rows: agg[(d['M'],d['K'])].append(d)
Ms=sorted(set(M for M,K in agg))
print('| M | X | K* (free fraction 1/2) | lambda* = K* log M | H = log^2 M/loglog M | lambda*/H | lambda*/(log X log M/loglog M) | pair degree at K* |')
print('|---|---|---|---|---|---|---|---|')
for M in Ms:
    Ks=sorted(K for (m,K) in agg if m==M)
    ff=[(K, sum(d['free_frac'] for d in agg[(M,K)])/len(agg[(M,K)]), sum(d['deg'] for d in agg[(M,K)])/len(agg[(M,K)])) for K in Ks]
    ks=None
    for (K1,f1,g1),(K2,f2,g2) in zip(ff,ff[1:]):
        if f1<0.5<=f2:
            t=(0.5-f1)/(f2-f1); ks=K1+t*(K2-K1); deg=math.exp(math.log(g1)+t*(math.log(g2)-math.log(g1))); break
    L=math.log(M); H=L*L/math.log(L); X=2*M
    lam=ks*L
    print(f'| {M} | {X} | {ks:.2f} | {lam:.1f} | {H:.1f} | {lam/H:.2f} | {lam/(math.log(X)*L/math.log(L)):.2f} | {deg:.1f} |')
