P=2147483647
class Ech:
    """incremental row echelon basis mod P over columns = keys"""
    def __init__(self): self.piv={}  # pivot col -> row dict (normalized pivot 1)
    def add(self,r):
        r={k:v%P for k,v in r.items() if v%P}
        while r:
            k=max(r)
            if k in self.piv:
                c=r[k]; pr=self.piv[k]
                for kk,vv in pr.items():
                    nv=(r.get(kk,0)-c*vv)%P
                    if nv: r[kk]=nv
                    else: r.pop(kk,None)
            else:
                inv=pow(r[k],P-2,P)
                self.piv[k]={kk:vv*inv%P for kk,vv in r.items()}
                return True
        return False
    def rank(self): return len(self.piv)
