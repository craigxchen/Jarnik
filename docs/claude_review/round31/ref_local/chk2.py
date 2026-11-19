from fractions import Fraction as F
def c(m): return m**3-11*m**2+39*m-13
def D(m): return m**4+70*m**2-48*m+9
def A(m): return -(m+3)*c(m)
def B(m): return 2*D(m)
def v2(n):
    k=0
    while n%2==0: n//=2;k+=1
    return k,n
# 2-adic: is n (nonzero 2-adic integer known mod 2^N) certainly a nonsquare?
def nonsq2(n,N):
    n%=2**N
    if n==0: return False
    k,o=v2(n)
    if k%2: return True
    # odd part known mod 2^(N-k)
    if N-k>=3: return o%8!=1
    if N-k==2: return o%4==3
    return False
tab={}
for r in range(1,256,2):
    tab.setdefault(r%8,set()).add((A(r)%256,B(r)%256))
for r in sorted(tab): print("r mod8",r,tab[r])
# depends on m mod 256 only since integer polys; check all odd residues mod 2^12 for robustness
ok=all(nonsq2(A(r),12) or nonsq2(B(r),12) for r in range(1,4096,2))
print("odd m: A or B nonsquare (mod 2^12 check):",ok)
ok=all(nonsq2(A(r),8) or nonsq2(B(r),8) for r in range(1,256,2)); print("mod 256 suffices:",ok)
# even m in Z_2: D odd
print("even m D odd:",all(D(r)%2==1 for r in range(0,256,2)))
# v2(m)<0: m=w/2^k, D = (w^4+70w^2 4^k... )/2^{4k}; numerator odd
print("v2(m)<0:",all((w**4+70*w**2*4**k-48*w*8**k+9*16**k)%2==1 for w in range(1,64,2) for k in range(1,5)))
# 3-adic split control: -A, -B squares
def v3(n):
    k=0
    while n%3==0: n//=3;k+=1
    return k,n
def nonsq3(n,N):
    n%=3**N
    if n==0: return False
    k,o=v3(n)
    if k%2: return True
    if N-k>=1: return o%3==2
    return False
res={}
for r in range(3**6):
    a=nonsq3(-A(r),6); b=nonsq3(-B(r),6)
    if not(a or b): print("FAIL",r)
print("all m in Z_3: -A or -B nonsquare")
print("m=1,2 mod3 -B mod3:",{(-B(r))%3 for r in range(3**5) if r%3}, " m=0 mod9 -A mod9:",{(-A(r))%9 for r in range(0,729,9)}," m=3 mod 9:",{(-A(r))%9 for r in range(3,729,9)}," m=6 mod9 -B mod81:",{(-B(r))%81 for r in range(6,729,9)})
# also: in Z_3 with m=6 mod 9, is -A possibly square? irrelevant
# twisted case at 3 for comparison: A,B squares possible?
print("twisted at 3 (A,B both possibly square) residues mod 27:",sorted({r%27 for r in range(729) if not(nonsq3(A(r),6) or nonsq3(B(r),6))})[:20])
# split at 2 for comparison
print("split at 2 surviving odd residues mod 8:",sorted({r%8 for r in range(4096) if not(nonsq2(-A(r),12) or nonsq2(-B(r),12))}))
