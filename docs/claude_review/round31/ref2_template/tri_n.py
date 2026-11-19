# exact check at one specialization: which signed triples of selected near-real ratios close (product = 1)
from itertools import combinations, product
from spec import B_at, C1, C2, SEL1, SEL2
def gm(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def cj(a): return (a[0],-a[1])
for name,roots,sel,n in (('C1',C1,SEL1,2**16*3**3*5**12*13**5*37**3*41**3*53**3*509**3+123457),
                         ('C2',C2,SEL2,10**40+987654321)):
    B=[B_at(roots,v,n) for v in sel]
    closing=[]
    for i,j,k in combinations(range(len(sel)),3):
        for e in product((1,-1),repeat=3):
            num=(1,0);den=(1,0)
            for idx,s in zip((i,j,k),e):
                x=B[idx] if s==1 else cj(B[idx]); y=cj(x)
                num=gm(num,x); den=gm(den,y)
            if num==den: closing.append((i,j,k,e))
    print(name,'n=',n,'closing signed triples (index into selected list):',closing)
