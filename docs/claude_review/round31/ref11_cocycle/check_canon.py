# verify A_j - A_i = 2 i eta H_ij T_ij with T_ij = +-Im X_ij, 0<T^2<Q, and content((A_j-A_i)/2)=|T|, on random pairs
import random, sys
from common import *
P = int(sys.argv[1]); ps, gs, q, s, marks = setup(P)
print("P",P,"s",s,"q",q, "q bits", q.bit_length())
random.seed(5); bad = 0
for _ in range(20000):
    i, j = random.sample(range(1<<s), 2)
    H = (1,0); X = (1,0)
    for k in range(s):
        a,b = gs[k]; fi = (a,-b) if (i>>k)&1 else (a,b)
        if ((i^j)>>k)&1: X = gmul(X, fi)
        else: H = gmul(H, fi)
    Q = Qof(ps,i,j); T = X[1]
    D = (marks[j][0]-marks[i][0], marks[j][1]-marks[i][1])
    # D should equal H * (conj X - X) = H * (-2 i Im X)
    R = gmul(H, (0, -2*T))
    ok = (D == R) and 0 < T*T < Q and Tof(marks[i], marks[j]) == abs(T) and X[0]**2+X[1]**2 == Q
    bad += not ok
print("bad", bad)
