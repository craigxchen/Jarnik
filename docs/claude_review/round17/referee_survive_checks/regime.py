from mpmath import mp, mpf, log, sqrt, ceil, floor, binomial, exp
mp.dps=50
import sympy
def params(q):
    d=q; r=2*q; n=q+1
    w=sqrt(8*d*log(4*r)); j=int(floor(w/2)); mprime=int(ceil(4*w))
    return d,r,n,w,j,mprime
q=3251
print('isprime',sympy.isprime(q),'q mod 4',q%4)
d,r,n,w,j,mp_=params(q)
print('316 ln 4r =',316*log(4*r),' d =',d, ' regime:', d>=316*log(4*r))
# smallest regime prime
for p in sympy.primerange(3000,3300):
    if p%4==3 and p>=316*log(8*p): print('first regime q',p); break
print('w',w,'j',j,"m'",mp_,"m'<q",mp_<q,"8j",8*j,"8j+8",8*j+8, "m'<=8j+8", mp_<=8*j+8)
def tauD(c0,C=1):
    return max((2/mpf(c0))*(2*log(4*r+5)+log(1000/mpf(C))), 4*log(4*r+5)+8, 2*log(1000*q*n**2/mpf(C)))
for c0 in [0.184,0.25,0.5,1,10]:
    t=tauD(c0)
    budget=q*t/2 - mp_*log(2*q) - log(q/mpf(mp_+1))
    print('c0',c0,'tau_D',t,'budget',budget,'budget/q',budget/q,'lnH\'=',mp_*q*t/20)
# tau at which budget = q ln 3 and = 0
for target,name in [(q*log(3),'q ln3'),(0,'0')]:
    t=2*(target+mp_*log(2*q)+log(q/mpf(mp_+1)))/q
    print('tau with budget=',name, t, ' (3/2)ln n =', 1.5*log(n), ' W/(n ln n) at that tau:', 2*q*t/(n*log(n)))
# W_min 3 n ln n - 14 n -> tau
Wmin=3*n*log(n)-14*n; print('W_min',Wmin,'tau=',Wmin/(2*q),'budget there',q*(Wmin/(2*q))/2-mp_*log(2*q)-log(q/mpf(mp_+1)))
# Hamming ball lower bound r0
s=sum(binomial(q,t)*2**t for t in range(j+1)); print('ln r0 >=',log(s))
