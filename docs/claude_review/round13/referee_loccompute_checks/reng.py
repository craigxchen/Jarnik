# Reduced ring R = Q[x_1..x_d]/(x^3), dense vectors mod p per degree; cyclic orbit projection P0.
import numpy as np, itertools
from rlib import monos, J_exact
PR = 2147483629
class REng:
    def __init__(self, d, Dmax, shift=None):
        self.d=d; self.Dmax=Dmax
        self.ms=[monos(d,D,2) for D in range(Dmax+1)]
        self.idx=[{m:i for i,m in enumerate(ms)} for ms in self.ms]
        # multiplication by x_j maps
        self.mul=[]
        for D in range(Dmax):
            arr=np.full((d,len(self.ms[D])),-1,dtype=np.int64)
            for i,m in enumerate(self.ms[D]):
                for j in range(d):
                    if m[j]<2:
                        mm=list(m); mm[j]+=1
                        arr[j,i]=self.idx[D+1][tuple(mm)]
            self.mul.append(arr)
        if shift is not None:  # shift: permutation of coordinates (list), orbit ids per degree
            self.orb=[]
            for D in range(Dmax+1):
                oid=np.full(len(self.ms[D]),-1,dtype=np.int64); cnt=0
                for i,m in enumerate(self.ms[D]):
                    if oid[i]>=0: continue
                    cur=m
                    while True:
                        k=self.idx[D][cur]
                        if oid[k]>=0: break
                        oid[k]=cnt
                        cur=tuple(cur[shift[j]] for j in range(d))
                    cnt+=1
                self.orb.append((oid,cnt))
    def one(self):
        return np.array([1],dtype=np.int64)
    def mul_lin(self, v, D, l):
        """v in R_D (mod p), l integer linear form (len d) -> R_{D+1}"""
        out=np.zeros(len(self.ms[D+1]),dtype=np.int64)
        for j in range(self.d):
            c=int(l[j])%PR
            if c==0: continue
            tgt=self.mul[D][j]; ok=tgt>=0
            np.add.at(out, tgt[ok], (v[ok]*c)%PR)
            out%=PR
        return out
    def mul_var(self, v, D, j):
        out=np.zeros(len(self.ms[D+1]),dtype=np.int64)
        tgt=self.mul[D][j]; ok=tgt>=0
        out[tgt[ok]]=v[ok]
        return out
    def pullback(self, coeffs, tmonos, U):
        """form sum c_e prod_i <U_i,x>^{e_i} (homogeneous of degree D) -> vector in R_D mod p"""
        cache={}
        def tv(e):
            if e in cache: return cache[e]
            if sum(e)==0: r=self.one()
            else:
                i=next(k for k in range(len(e)) if e[k]>0)
                ee=list(e); ee[i]-=1; ee=tuple(ee)
                r=self.mul_lin(tv(ee), sum(ee), U[i])
            cache[e]=r; return r
        D=sum(tmonos[0]); out=np.zeros(len(self.ms[D]),dtype=np.int64)
        for c,e in zip(coeffs,tmonos):
            if c: out=(out+ (int(c)%PR)*tv(e))%PR
        return out
    def ideal_rows(self, g, D, Dtarget):
        """all mu*g, mu monomial of R_{Dtarget-D}; yields vectors in R_Dtarget"""
        # BFS over monomials via variable multiplication, avoiding duplicates by nondecreasing var index
        res=[]
        def rec(v,DD,start):
            if DD==Dtarget: res.append(v); return
            for j in range(start,self.d):
                rec(self.mul_var(v,DD,j),DD+1,j)
        rec(g,D,0)
        return res
    def P0(self, v, D):
        oid,cnt=self.orb[D]
        out=np.zeros(cnt,dtype=np.int64)
        np.add.at(out,oid,v); return out%PR
