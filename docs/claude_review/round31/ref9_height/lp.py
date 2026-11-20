from scipy.optimize import linprog
# vars xi,ka,ta,a ; maximize xi; closed versions of the strict constraints
A=[[2,-1,0,0],[1,1,-1,0],[0,0,1,0],[1,1,0,1]]; bb=[0,0,0.5,0.5]
r=linprog(c=[-1,0,0,0],A_ub=A,b_ub=bb,bounds=[(0,None)]*4); print('new system sup xi =',-r.fun,r.x)
r=linprog(c=[-1,0,0,0],A_ub=A+[[0,0,1,-1]],b_ub=bb+[0],bounds=[(0,None)]*4); print('old system (a>=tau) sup xi =',-r.fun,r.x)
