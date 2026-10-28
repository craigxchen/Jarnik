import mpmath as mp
mp.mp.dps=40
from c4_construct import curve, funcs
def angles(m,w):
    C=curve(m); Xz,P,H,U,kap,sk=funcs(C)
    return [mp.arg(x) for x in U(w)]
for m0,w0 in ((-1.5,0.7),(-2.2,-1.3),(-2.7,0.05),(-1.2,3.0)):
    m0=mp.mpf(m0); w0=mp.mpf(w0); h=mp.mpf(10)**-15
    a=angles(m0,w0)
    dm=[(x-y)/(2*h) for x,y in zip(angles(m0+h,w0),angles(m0-h,w0))]
    dw=[(x-y)/(2*h) for x,y in zip(angles(m0,w0+h),angles(m0,w0-h))]
    J=mp.matrix([dm,dw]).T
    print(float(m0),float(w0),'singular values',[mp.nstr(s,4) for s in mp.svd_r(J,compute_uv=False)])
