# Witness-height budget of Corollary 3.2 at regime sdc_q (C=1). 50-digit mpmath.
from mpmath import mp, mpf, sqrt, log, ceil, floor, binomial
mp.dps=50
def isprime(n):
    if n<2: return False
    i=2
    while i*i<=n:
        if n%i==0: return False
        i+=1
    return True
def params(q):
    d=q; r=2*q; n=q+1
    w=sqrt(8*d*log(4*r)); j=int(floor(w/2)); mp_=int(ceil(4*w))
    return d,r,n,w,j,mp_
def tauD(q,c0,C=1):
    d,r,n,w,j,m=params(q)
    return max((2/mpf(c0))*(2*log(4*r+5)+log(1000/mpf(C))), 4*log(4*r+5)+8, 2*log(1000*q*n**2/mpf(C)))
def budget(q,tau,C=1):
    d,r,n,w,j,m=params(q)
    return q*tau/2 - m*log(2*q) - log(C*mpf(q)/(m+1))
if __name__=="__main__":
  qs=[q for q in range(3251,3400) if isprime(q) and q%4==3][:3]+[5479,10007,100003]
  qs=[q for q in qs if isprime(q) and q%4==3]
  print("q, j, m', m'<q, tau_D(c0=.184), budget lnL*(tau_D), budget/q, q ln3, tau_min for budget>=q ln 3, tau_min for budget>=0")
  for q in qs:
      d,r,n,w,j,m=params(q)
      tD=tauD(q,0.184)
      B=budget(q,tD)
      # tau with budget = q ln 3
      t3=2*(q*log(3)+m*log(2*q)+log(mpf(q)/(m+1)))/q
      t0=2*(m*log(2*q)+log(mpf(q)/(m+1)))/q
      print(q,j,m,m<q,float(tD),float(B),float(B/q),float(q*log(3)),float(t3),float(t0))
  for c0 in [0.184,0.1,0.05]:
      print("c0",c0,"tau_D(3251)=",float(tauD(3251,c0)),"budget/q=",float(budget(3251,tauD(3251,c0))/3251))
