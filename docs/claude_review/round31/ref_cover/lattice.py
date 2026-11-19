def lattice_index(cols):
    """cols: list of integer 4-vectors. Return (rank, index or None) of the Z-span in Z^4, by integer row echelon (HNF-like)."""
    rows=[list(c) for c in cols if any(c)]
    basis=[]  # echelon basis
    n=4
    piv_rows=[]
    M=rows
    rank=0; det=1
    col=0
    while col<n and M:
        nz=[r for r in M if r[col]!=0]
        z=[r for r in M if r[col]==0]
        if not nz:
            col+=1; continue
        # gcd-reduce among nz on column col
        while len(nz)>1:
            nz.sort(key=lambda r: abs(r[col]))
            p=nz[0]; new=[p]
            for r in nz[1:]:
                q=r[col]//p[col]
                rr=[a-q*b for a,b in zip(r,p)]
                if rr[col]!=0: new.append(rr)
                elif any(rr): z.append(rr)
            nz=new
        p=nz[0]
        det*=abs(p[col]); rank+=1
        M=z; col+=1
    if rank<4: return rank,None
    return 4,det
if __name__=="__main__":
    import random, sympy as sp
    from itertools import combinations
    from math import gcd
    for t in range(2000):
        k=random.randint(1,7)
        cols=[[random.randint(-3,3) for _ in range(4)] for _ in range(k)]
        r,I=lattice_index(cols)
        M=sp.Matrix(cols).T
        assert r==M.rank()
        if r==4:
            g=0
            for c in combinations(range(k),4): g=gcd(g,int(abs(M[:,list(c)].det())))
            assert g==I,(cols,g,I)
    print("lattice_index ok")
