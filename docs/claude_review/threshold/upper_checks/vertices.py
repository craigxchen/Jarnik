from fractions import Fraction as F
def T(x): return x*(x+1)//2
def cons(M):
    return [(a,b) for a in range(0,M) for b in range(0,M) if a+b<=M-1 and a+b>=1]
def feasible(M, al, be):
    return all(al*T(b)+be*T(a) >= T(a+b) for (a,b) in cons(M))
def vertices(M):
    C = cons(M); V=set()
    for i in range(len(C)):
        for j in range(i+1,len(C)):
            a1,b1=C[i]; a2,b2=C[j]
            # al*T(b1)+be*T(a1)=T(a1+b1); al*T(b2)+be*T(a2)=T(a2+b2)
            det = T(b1)*T(a2)-T(b2)*T(a1)
            if det==0: continue
            al = F(T(a1+b1)*T(a2)-T(a2+b2)*T(a1), det)
            be = F(T(b1)*T(a2+b2)-T(b2)*T(a1+b1), det)
            if al>=1 and be>=1 and feasible(M,al,be):
                tight = [(a,b) for (a,b) in C if al*T(b)+be*T(a)==T(a+b)]
                V.add((al,be,tuple(tight)))
    return sorted(V)
for M in range(2,10):
    print("M=",M)
    for al,be,tight in vertices(M):
        print("   ",al,be," tight (a,b):",tight)
