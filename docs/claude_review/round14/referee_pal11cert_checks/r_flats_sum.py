from collections import defaultdict
mint={}; att=defaultdict(set); dims=defaultdict(set); lines=defaultdict(set)
for ln in open('r_flats_out.txt'):
    a=list(map(int,ln.split())); R,Z,s,t,dn=a[:5]
    mint[s]=min(mint.get(s,99),t); att[s].add(t); dims[(s,t)].add(dn)
    if dn==1:
        v=a[5:]
        # normalise: divide by first nonzero, keep as fraction mod p? entries are small ints mod p after scaling by gcd? use sign/gcd
        from math import gcd
        g=0
        for x in v: g=gcd(g,abs(x))
        v=[x//g for x in v]
        nz=[x for x in v if x][0]
        if nz<0: v=[-x for x in v]
        lines[(s,t)].add(tuple(v))
print('min t by sigma', [ (s,mint[s]) for s in sorted(mint)])
for s in sorted(att): print('sigma',s,'attainable t',sorted(att[s]))
for k in sorted(dims):
    if k[1]<=11: print(k,'null dims',sorted(dims[k]), 'nlines(dim1)',len(lines[k]), 'maxentry', max((max(abs(x) for x in v) for v in lines[k]),default=None))
