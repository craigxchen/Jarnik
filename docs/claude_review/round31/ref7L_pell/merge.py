import itertools
F={'G1':[0,1,1,0,1],'G1b':[1,0,0,1,0],'G2':[0,1,0,1,0],'G2b':[1,0,1,0,1],'G3':[0,0,1,1,0],'G3b':[1,1,0,0,1],
   'G4':[0,0,0,0,1],'G4b':[1,1,1,1,0],'b':[1,0,0,0,0],'bb':[0,1,1,1,1]}
for (p,P),(q,Q) in itertools.combinations(F.items(),2):
    S=[a+b for a,b in zip(P,Q)]; red=(max(P)-min(P))+(max(Q)-min(Q))-(max(S)-min(S))
    if red>0: print(p,q,S,"reduction per shared root",red)
