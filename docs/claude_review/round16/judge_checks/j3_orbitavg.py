# Judge check J3 (50-digit numerics; the statement itself has a 2-line proof): for real f homogeneous of degree m
# in (cos t, sin t), the average of f^2 over N equispaced points (offset t0, step 2*pi*a/N, gcd(a,N)=1) equals
# the mean of f^2 over the circle whenever N is odd and m < N.  Fails in general at m = N.
import mpmath as mp, random
mp.mp.dps=50; random.seed(1)
def test(m,N,a=1):
    co=[mp.mpf(random.randint(-9,9)) for _ in range(m+1)]
    f=lambda t: sum(co[i]*mp.cos(t)**i*mp.sin(t)**(m-i) for i in range(m+1))
    mean=mp.quad(lambda t: f(t)**2,[0,mp.pi/2,mp.pi,3*mp.pi/2,2*mp.pi])/(2*mp.pi)
    t0=mp.mpf('0.3217')
    avg=sum(f(t0+2*mp.pi*a*s/N)**2 for s in range(N))/N
    return avg/mean
for N in (7,11):
    for m in range(1,N+2):
        print("N=%d m=%d ratio=%s"%(N,m,mp.nstr(test(m,N,a=2),20)))
