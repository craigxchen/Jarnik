from spec import *
seedsets=[(F(1,2),F(-2),F(7,3),F(4)),(F(-1,3),F(5,2),F(6),F(-7)),(F(0),F(1),F(3,5),F(-11,4))]
for n in (1,2,3):
  for sd in seedsets:
    for T in range(-50,200):
        zs=tuple_at(n,sd,T)
        if len(set(zs))<5:
            print(n,[str(x) for x in sd],T,[ (a,b) for a in range(5) for b in range(a+1,5) if zs[a]==zs[b]])
