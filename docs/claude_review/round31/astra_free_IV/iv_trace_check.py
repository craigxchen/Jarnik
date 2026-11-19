"""Exact positive-characteristic control for the polynomial integral-fake model."""
from itertools import combinations


def check(degree, modulus):
    q = 1 << degree

    def mul(a,b):
        s = 0
        while b:
            if b & 1:
                s ^= a
            b >>= 1
            a <<= 1
            if a & q:
                a ^= modulus
        return s

    def power(a,n):
        s=1
        while n:
            if n & 1:
                s=mul(s,a)
            a=mul(a,a)
            n//=2
        return s

    assert all(power(a,q-1)==1 for a in range(1,q))

    def trace(a):
        s=0
        for _ in range(degree):
            s ^= a
            a = mul(a,a)
        assert s in (0,1)
        return s

    def multiply_linear(p,x):
        result=[0]*(len(p)+1)
        for i,c in enumerate(p):
            result[i] ^= mul(x,c)
            result[i+1] ^= c
        return result

    rows={}
    for a in range(1,q):
        marked={x for x in range(1,q) if trace(mul(a,x))==0}
        assert len(marked)==q//2-1
        coeff=[0]*(q//4)
        for j in range(degree-1):
            n=1<<j
            coeff[n-1]=power(a,(n-q//2)%(q-1))
        rows[a]=marked,coeff
    for i,j in combinations(rows,2):
        si,ci=rows[i]
        sj,cj=rows[j]
        assert len(si^sj)==q//2
        g=[1]
        for x in si&sj:
            g=multiply_linear(g,x)
        det=[a^b for a,b in zip(ci,cj)]
        assert det[-1]!=0
        assert len(g)==len(det) and all(v==mul(w,det[-1]) for v,w in zip(det,g))
    print('q, rows, pairs with nonzero constant cofactor, pair separation:',
          q,q-1,(q-1)*(q-2)//2,q//2)


if __name__=='__main__':
    for d,modulus in ((3,0b1011),(4,0b10011),(5,0b100101),(7,0b10000011)):
        check(d,modulus)
