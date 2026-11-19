# Independent implementation of the III.10 one-list test over F_p (split chart).
import sys, itertools
from sympy import Matrix, GF
from sympy.polys.matrices import DomainMatrix
def run(p, roots):
    edges=[(i,j) for i in range(5) for j in range(i+1,5)]
    u=dict(zip(edges,roots))
    F=GF(p)
    def pmul(A,B):
        C=[0]*(len(A)+len(B)-1)
        for i,a in enumerate(A):
            for j,b in enumerate(B): C[i+j]=(C[i+j]+a*b)%p
        return C
    P=[]
    for i in range(5):
        f=[1]
        for e in edges:
            if i in e: f=pmul(f,[(-u[e])%p,1])   # low->high coefficients
        P.append(f)
    M=DomainMatrix([[F(P[i][k]) for k in range(5)] for i in range(5)],(5,5),F)
    rk=M.rank(); print("p",p,"rank",rk)
    ns=M.nullspace().to_Matrix()
    if rk!=4: print("not rank 4"); return
    a=[int(x)%p for x in ns.row(0)]  # a_1..a_5
    # R_i = polynomial part of a(x) P_i(x), a(x)=sum a_k x^{-k}
    R=[]
    for i in range(5):
        r=[0]*4
        for k in range(1,6):
            for m in range(5):
                if m-k>=0: r[m-k]=(r[m-k]+a[k-1]*P[i][m])%p
        R.append(r)
        assert sum(P[i][m]*a[m] for m in range(5))%p==0
    def ev(f,x): return sum(c*pow(x,k,p) for k,c in enumerate(f))%p
    lam={};v={}
    for (i,j) in edges:
        num=[(x-y)%p for x,y in itertools.zip_longest(pmul(P[i],R[j]),pmul(P[j],R[i]),fillvalue=0)]
        # divide by (x-u_ij)
        while num and num[-1]==0: num.pop()
        # synthetic division
        r=u[(i,j)]; q=[0]*(len(num)-1); acc=0
        for k in range(len(num)-1,-1,-1):
            acc=(acc*r+num[k])%p
            if k>0: q[k-1]=acc
        assert acc==0
        while q and q[-1]==0: q.pop()
        assert len(q)==2, q
        vij=(-q[0]*pow(q[1],-1,p))%p; v[(i,j)]=vij
        li=(-ev(P[i],vij)*pow(ev(R[i],vij),-1,p))%p
        lj=(-ev(P[j],vij)*pow(ev(R[j],vij),-1,p))%p
        assert li==lj
        lam[(i,j)]=li
    allr=list(u.values())+list(v.values())
    print("twenty distinct:",len(set(allr))==20)
    print("v:",[v[e] for e in edges]); print("lambda:",[lam[e] for e in edges])
run(101,[35,76,91,75,13,88,93,5,3,92])
