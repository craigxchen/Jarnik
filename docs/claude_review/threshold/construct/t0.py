from clib import *
from gens import *
for n in range(1,5):
    S = Qslice(3,n,n); print("Q3",n,n,len(S),reg2(proj(S)))
for M in (3,4,5):
    r=(M-1)//2
    v=tuple(range(-r, M-r))
    S=orbit(v); print("orbit",v,len(S),reg2(proj(S)), "C(M,2)=",M*(M-1)//2)
print("Q4(3,1)",reg2(proj(Qslice(4,3,1))))
print("Q5(3,2)",reg2(proj(Qslice(5,3,2))), "Q5(3,3)", reg2(proj(Qslice(5,3,3))))
# LP test: M=3 widths
w={1:4,2:2,3:3}
print(adversary(w,3))
