import sys, itertools
def deg(a,b,c):
    m=max(a,b,c); return 2*m==a+b+c
def line_ok(n,d):
    # find diametral pair, coordinates from one end
    u,v=max(((i,j) for i in range(n) for j in range(i+1,n)),key=lambda t:d[t[0]][t[1]])
    x=[d[u][i] for i in range(n)]
    return all(abs(x[i]-x[j])==d[i][j] for i in range(n) for j in range(n))
def pseudolinear(d):
    # 4 points: some labelling with d01=d23, d12=d30, d02=d13=d01+d12
    for p in itertools.permutations(range(4)):
        a,b,c,e=p
        if d[a][b]==d[c][e] and d[b][c]==d[e][a] and d[a][c]==d[b][e]==d[a][b]+d[b][c]: return True
    return False
def run(n,N):
    d=[[0]*n for _ in range(n)]
    cnt=dict(total=0,line=0,pseudo=0,other=0)
    others=[]
    # fix vertex order: add vertex k with distances to 0..k-1
    def rec(k):
        if k==n:
            cnt['total']+=1
            if line_ok(n,d): cnt['line']+=1
            elif n==4 and pseudolinear(d): cnt['pseudo']+=1
            else:
                cnt['other']+=1
                if len(others)<5: others.append([r[:] for r in d])
            return
        def rec2(j):
            if j==k: rec(k+1); return
            for t in range(1,N+1):
                d[k][j]=d[j][k]=t
                if all(deg(d[i][j],d[i][k],t) for i in range(j)):
                    rec2(j+1)
            d[k][j]=d[j][k]=0
        rec2(0)
    rec(1)
    print(n,N,cnt,others[:3])
for n,N in [(4,24),(5,20),(6,14)]:
    run(n,N)
