from mpmath import mp, mpf, atan2, sqrt, pi, sin, cos, fabs, floor, atan
mp.dps=60
U=pi/2
def red(a):  # distance to (pi/2)Z
    a=a % U
    return min(a, U-a)
# three fibres p=29
phi=atan(sqrt(2))
beta=[(533,8422),(1613,994),(139,326)]
hs=[29**j*(a*a+b*b) for j,(a,b) in enumerate(beta)]
print('heights',hs)
X=max(hs)
angs=[(atan2(b,a)+j*phi)%U for j,(a,b) in enumerate(beta)]
from itertools import combinations
w=max(red(a-b) for a,b in combinations(angs,2))
print('X*width =', mp.nstr(X*w,10))
for a,b in beta:
    from math import gcd; assert gcd(a,b)==1 and (a*a+b*b)%29!=0
# monomial gaps for phi=arctan sqrt2 at p=29 : n up to 200
worst_ax=min(sin(red(n*phi))*mpf(29)**(mpf(n)/2) for n in range(1,201))
worst_dg=min(sin(red(n*phi-pi/4))*mpf(29)**(mpf(n)/2) for n in range(1,201))
print('p=29: min_n sin||n phi||*29^{n/2} =',mp.nstr(worst_ax,6),' diagonal:',mp.nstr(worst_dg,6), '(claimed >=1, >=1/sqrt2)')
# p=13 graft phi=theta/2
th=atan2(3,2); ph=th/2
print('tan(phi) vs (sqrt13-2)/3:', mp.nstr(mp.tan(ph)-(sqrt(13)-2)/3,5))
ax=min(sin(red(n*ph))*mpf(13)**(mpf(n)/2) for n in range(1,301))
dg=min(sin(red(n*ph-pi/4))*mpf(13)**(mpf(n)/2) for n in range(1,301))
print('p=13: min sin||n phi||13^{n/2}=',mp.nstr(ax,6),' diag:',mp.nstr(dg,6))
