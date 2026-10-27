# necessary condition K*o_q > d for some admissible q (bit strata), n=2^k, r=n-1
import mpmath as mp
mp.mp.dps=40
def cmin(q,n,om):
    # least integer c>=0 with q*2*exp(-3(c+1)^2/(2n)) <= om ; solve then adjust
    c=max(0,int(mp.sqrt(2*n*mp.log(2*q/om)/3))-3)
    while q*2*mp.e**(-3*mp.mpf(c+1)**2/(2*n))>om: c+=1
    while c>0 and q*2*mp.e**(-3*mp.mpf(c)**2/(2*n))<=om: c-=1
    return c
res={}
for th in ['half','t0','tS']:
    last=None
    for k in range(12,61):
        n=2**k; r=n-1; d=n-1; w=mp.sqrt(8*d*mp.log(4*r)); m=int(mp.ceil(4*w))
        om={'half':mp.mpf(1)/2,'t0':mp.mpf(1)/2*(4*r)**mp.mpf(-0.5),'tS':2*r*mp.e**(-3*w**2/(8*n))}[th]
        ok=False; qmax=0
        for q in range(1,k+1):
            c=cmin(q,n,om); rho=2*c*q
            if rho<=m-1:
                qmax=q
                if 2**q*(m-rho)>d: ok=True
        if ok: last=k
        if th=='tS' and k in (12,30,60): print('tS k',k,'qmax',qmax)
    print(th,'largest k in [12,60] with K*o_q>d for some q:',last)
print('detail tS:')
for k in (12,13,14):
    n=2**k; r=n-1; d=n-1; w=mp.sqrt(8*d*mp.log(4*r)); m=int(mp.ceil(4*w)); j=int(mp.floor(w/2))
    om=2*r*mp.e**(-3*w**2/(8*n))
    for q in range(1,6):
        c=cmin(q,n,om); rho=2*c*q
        print(k,'q',q,'c',c,'j',j,'rho',rho,'m',m,'o',m-rho,'K*o',2**q*(m-rho),'d',d, 'adm',rho<=m-1)
