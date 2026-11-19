# exponent tables and common degree g for q=3 (labels = F_3^k) and q=27,k=2
import itertools, sys
def chi3(a): a%=3; return 0 if a==0 else (1 if a==1 else -1)
def forms_proj(k,p=3):
    out=[]
    for c in itertools.product(range(p),repeat=k):
        if any(c):
            f=[x for x in c if x][0]
            if f==1: out.append(c)
    return out
def table(k, forms, signs):
    rows=[]
    labs=list(itertools.product(range(3),repeat=k))
    for f,e in zip(forms,signs):
        rows.append([1+chi3(e*sum(a*b for a,b in zip(f,x))) for x in labs])
    return labs, rows
def g_of(k,forms,signs):
    labs,rows=table(k,forms,signs)
    return sum(min(c) for c in zip(*rows))
for k in (2,3):
    n=3**k; h=n//3; P=forms_proj(k)
    print("k",k,"n",n,"h",h,"#forms",len(P))
    for m in range(2,len(P)+1):
        best=None
        for fs in itertools.combinations(P,m):
            for sg in itertools.product((1,-1),repeat=m):
                if sg[0]==-1: continue
                g=g_of(k,fs,sg)
                if best is None or g>best[0]: best=(g,fs,sg)
            if m>6: break
        print(" m",m,"best g",best[0],"g-h",best[0]-h,best[1],best[2])
