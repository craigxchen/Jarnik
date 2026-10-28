# Exhaustive sanity check of (R_n), n=3,4,5 (Lemma 2.4) and of monotonicity of F_n, for 0<=a,b<=L.
import sys
L=int(sys.argv[1])
H={2:lambda x,y:0,3:lambda x,y:min(x,y),4:lambda x,y:min((x+y)//2,2*min(x,y)),
   5:lambda x,y:min(min(x,y)+max(x,y)//3,3*min(x,y))}
for n in (3,4,5):
    ok=all(min(H[n-1](a-p,b)+q,H[n-1](a,b-q)+p)<=H[n](a,b) for a in range(L+1) for b in range(L+1) for p in range(a+1) for q in range(b+1))
    # tightness: is H_n the max of the left side?
    tight=all(max(min(H[n-1](a-p,b)+q,H[n-1](a,b-q)+p) for p in range(a+1) for q in range(b+1))==H[n](a,b) for a in range(L+1) for b in range(L+1))
    mono=all(a+1+b+H[n](a+1,b)>=a+b+H[n](a,b) and a+b+1+H[n](a,b+1)>=a+b+H[n](a,b) for a in range(L) for b in range(L))
    print(f"n={n}: (R_n) holds on [0,{L}]^2: {ok}; H_n equals max_(p,q) of left side: {tight}; F_n monotone: {mono}")
