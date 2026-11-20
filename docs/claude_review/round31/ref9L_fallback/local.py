# Exhaustive local check: for p in {2,3,5,7}, 7 finite marks, at most one T (3<=|T|<=6) with p | e_T,
# v_p(e_T)=a; marks of T congruent mod p^a. Count pairs i<j with p | c_ij,
# c_ij = (y_i-y_j)/prod_{T'>= {i,j}} e_T'.  Only the p-part matters.
# Minimize over configurations.  Also: count distinct visible splits T (2<=|T|<=6) with l_T(p)>0
# for random primitive 7-tuples (contact integers of the actual point).
import itertools, random, sys
def vp(n,p):
    if n==0: return 10**9
    v=0
    while n%p==0: n//=p; v+=1
    return v
def mincount(p, a_max=2):
    best=(10**9,None)
    subsets=[None]+[T for r in range(3,7) for T in itertools.combinations(range(7),r)]
    for T in subsets:
        for a in ([0] if T is None else range(1,a_max+1)):
            mod=p**(a+1)
            # marks in T: fix T's residue class r=0 mod p^a WLOG (translation); lifts mod p^(a+1): p choices
            # marks outside T: residues mod p^(a+1)
            Tset=set(T) if T else set()
            choices=[]
            for i in range(7):
                if i in Tset: choices.append([k*p**a for k in range(p)])
                else: choices.append(list(range(mod)))
            # reduce: symmetric search by enumerating (too big for p=7,a=2); restrict a<=1 for p>=5
            cnt_best=10**9
            for ys in itertools.product(*choices):
                c=0
                for i,j in itertools.combinations(range(7),2):
                    d=(ys[i]-ys[j])%mod
                    need = a if (i in Tset and j in Tset) else 0
                    # p | c_ij iff v_p(y_i-y_j) > need; with ys mod p^(a+1) this is d % p^(need+1)==0
                    if d % (p**(need+1))==0: c+=1
                    if c>=cnt_best: break
                if c<cnt_best: cnt_best=c
                if cnt_best==0: break
            if cnt_best<best[0]: best=(cnt_best,(T,a))
    return best
if __name__=="__main__":
    for p,amax in [(2,3),(3,2)]:
        print("p",p,"min #cofactors divisible by p:",mincount(p,amax)); sys.stdout.flush()
