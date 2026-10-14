from rlib import *
import pickle, random, itertools
A,H,dimJ,W=pickle.load(open('m8_a.pkl','rb'))
p=P; random.seed(2)
def image(V): return sorted(set(tuple(sum(a*b for a,b in zip(k,v)) for v in V) for k in A))
lines=[(a,b,a^b) for a in range(1,8) for b in range(a+1,8) if (a^b)>b]
print("Fano lines",lines)
fam={}
for F in itertools.combinations(range(1,8),2): fam[('pair',F)]=[col(a) for a in F]
for F in itertools.combinations(range(1,8),3):
    typ='line' if (F[0]^F[1])==F[2] else 'nonline'
    fam[(typ,F)]=[col(a) for a in F]
info={}
done={}
for key,V in fam.items():
    Y=image(V); b=len(V)
    sig=(key[0],len(Y))
    r=rho(Y,b)
    g=leading_ideal_gens(Y,b,8)
    info[key]=(len(Y),r,{e:len(g[e]) for e in g})
    done[key]=g
    if sig not in [s for s in [(k[0],info[k][0]) for k in list(info)[:-1]]]:
        print(key,"|Y|",len(Y),"rho",r,"dims J(Y)_e",[len(g[e]) for e in range(9)],flush=True)
        if key[0]=='line':
            full=[(y1,y2,y3) for y1 in range(-2,3) for y2 in range(-2,3) for y3 in range(-2,3) if (y1+y2+y3)%2==0]
            print("   parity box size",len(full),"missing",sorted(set(full)-set(Y)), "Y subset parity:",set(Y)<=set(full))
print("summary types:",set((k[0],v[0],v[1]) for k,v in info.items()))
pickle.dump((fam,done),open('m8_c.pkl','wb'))
