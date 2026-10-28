# Independent check of Lemmas 2.3-2.4 and the closed recursion (referee code; not reusing width_checks).
import sys
def H(n,x,y):
    lo,hi=min(x,y),max(x,y)
    if n==2: return 0
    if n==3: return lo
    if n==4: return min((x+y)//2, 2*lo)
    if n==5: return min(lo+hi//3, 3*lo)
def F(n,x,y): return x+y+H(n,x,y)
