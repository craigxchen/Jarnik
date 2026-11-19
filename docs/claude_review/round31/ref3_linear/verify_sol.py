# direct check of a printed control solution: Im power sums of the five stars, k=1..5
import re,sys
def chk(p,line):
    stars=[ [tuple(map(int,re.match(r'(\d+)\+(\d+)i',t).groups())) for t in s.split()] for s in re.findall(r'\{([^}]*)\}',line)]
    def mul(a,b): return ((a[0]*b[0]-a[1]*b[1])%p,(a[0]*b[1]+a[1]*b[0])%p)
    def pw(a,k):
        r=(1,0)
        for _ in range(k): r=mul(r,a)
        return r
    allel=[e for s in stars for e in s]; from collections import Counter
    c=Counter(allel); assert all(v==2 for v in c.values()) and len(c)==10
    canon={min(e,(e[0],(-e[1])%p)) for e in c}; assert len(canon)==10 and all(e[1]!=0 for e in c)
    for k in range(1,6):
        print(p,"k",k,[sum(pw(e,k)[1] for e in s)%p for s in stars])
chk(7,"{0+1i 2+1i 4+6i 5+6i } {0+1i 1+3i 0+5i 5+5i } {2+1i 1+3i 2+5i 4+5i } {0+5i 1+5i 4+5i 4+6i } {1+5i 2+5i 5+5i 5+6i }")
chk(11,"{1+1i 7+3i 2+4i 9+4i } {1+1i 5+4i 0+8i 10+10i } {10+2i 2+4i 10+7i 10+10i } {10+2i 9+4i 0+8i 6+9i } {7+3i 5+4i 10+7i 6+9i }")
