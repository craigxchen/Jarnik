from rlib import *
import pickle, time
A,H,dimJ,W=pickle.load(open('m8_a.pkl','rb'))
deg=8; basis=monos(7,deg); idx={m:i for i,m in enumerate(basis)}
rows=[]
for v,t,h in W:
    if h+1>deg: continue
    f=ppow(linform(v),h+1)
    for m in monos(7,deg-h-1): rows.append(vec(pmul(f,{m:1}),idx))
print("rows",len(rows),flush=True)
t0=time.time(); r=fmpz_mat(rows).rank(); print("deg 8 EXACT Q rank (J_w)_8:",r,"time",round(time.time()-t0),flush=True)
