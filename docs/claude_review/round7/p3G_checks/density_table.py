"""Summarise density_out.txt (units s = 1, 2, 3; s = 0 has D_0 even, never tight, Remark 1.5)."""
import re
cur=None; rows={}
for line in open('density_out.txt'):
    m=re.match(r'M=(\d+) K=(\d+) X=(\d+)',line)
    if m: cur=(int(m.group(1)),int(m.group(2))); rows[cur]={'S':[0,0,0.0],'P':[0,0,0.0]}; continue
    m=re.match(r"\s+\('(S|P)', (.*)\) (\d+) (\d+) ([\d.]+)",line)
    if m and cur:
        kind=m.group(1); s=int(m.group(2).split(',')[-1])
        if s==0: continue
        n,o,p=int(m.group(3)),int(m.group(4)),float(m.group(5))
        rows[cur][kind][0]+=n; rows[cur][kind][1]+=o; rows[cur][kind][2]+=p
print('| M | K | X | singles: observed / predicted | rate | pairs: observed / predicted | rate |')
print('|---|---|---|---|---|---|---|')
for (M,K),v in sorted(rows.items()):
    S=v['S'];P=v['P']
    print(f"| {M} | {K} | {2*M} | {S[1]} / {S[2]:.0f} | {S[1]/S[0]:.3f} | {P[1]} / {P[2]:.0f} | {P[1]/P[0]:.4f} |")
