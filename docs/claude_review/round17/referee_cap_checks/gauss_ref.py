# Referee: exact check of (1-g) khat(1) = 2 sum_a zeta^a S_a, (1+g) khat(1) = 2 sum_a zeta^a T_a in Z[x]/Phi_q,
# with patterns built from the record definition (sdc_sym.py style): P_a={a}u{x: chi(a-x)=1}, P'_a={a}u{x: chi(a-x)=-1}.
import flint, random
def chi(a,q):
    a%=q
    return 0 if a==0 else (1 if pow(a,(q-1)//2,q)==1 else -1)
def primes3():
    for q in range(7,140):
        if q%4==3 and all(q%p for p in range(2,q)): yield q
random.seed(17)
allok=True
for q in primes3():
    Phi=flint.fmpz_poly([1]*q)
    X=lambda e: flint.fmpz_poly([0]*(e%q)+[1])
    red=lambda f: f % Phi
    g=red(sum((X(x)*chi(x,q) for x in range(1,q)), flint.fmpz_poly([0])))
    ok_g = red(g*g)==flint.fmpz_poly([-q])
    P=[set([a]+[x for x in range(q) if chi(a-x,q)==1]) for a in range(q)]   # record P_a = {a} u (a+NQR)
    Pp=[set([a]+[x for x in range(q) if chi(a-x,q)==-1]) for a in range(q)] # record P'_a = {a} u (a+QR)
    nqr={x for x in range(1,q) if chi(x,q)==-1}; qr={x for x in range(1,q) if chi(x,q)==1}
    ok_def = all(P[a]=={a}|{(a+y)%q for y in nqr} and Pp[a]=={a}|{(a+y)%q for y in qr} for a in range(q))
    ok=True
    for trial in range(25):
        k=[random.randint(-5,5) for _ in range(q)]
        kh=red(sum((X(x)*k[x] for x in range(q)), flint.fmpz_poly([0])))
        S=[sum(k[x] for x in Pp[a]) for a in range(q)]
        T=[sum(k[x] for x in P[a]) for a in range(q)]
        lhs1=red((flint.fmpz_poly([1])-g)*kh); rhs1=red(2*sum((X(a)*S[a] for a in range(q)),flint.fmpz_poly([0])))
        lhs2=red((flint.fmpz_poly([1])+g)*kh); rhs2=red(2*sum((X(a)*T[a] for a in range(q)),flint.fmpz_poly([0])))
        K=sum(k)
        ok = ok and lhs1==rhs1 and lhs2==rhs2 and all(S[a]+T[a]==K+k[a] for a in range(q))
    print(f"q={q}: g^2=-q {ok_g}; patterns match author's P_a,P'_a {ok_def}; identities (2.1) and S+T=K+k_a on 25 random k: {ok}",flush=True)
    allok = allok and ok and ok_g and ok_def
print("ALL OK:",allok)
