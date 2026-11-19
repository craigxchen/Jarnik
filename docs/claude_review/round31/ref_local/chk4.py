exec(open('chk2.py').read().split('tab={}')[0])
from collections import Counter
cnt=Counter()
for x in range(5,2**14,8):
    a=nonsq2(-A(x),14); b=nonsq2(-B(x),14)
    cnt[(x%32, a, b)]+=1
for k in sorted(cnt): print(k,cnt[k])
for x in (5,13,21,29):
    print(x, v2(-A(x)), v2(-B(x)))
# twisted at 3 independent: for each m mod 3^6 report reason
def v3(n):
    k=0
    while n%3==0: n//=3;k+=1
    return k,n
bad=[x for x in range(3**7) if not((v3(A(x))[0]%2==1 and A(x)%3**7!=0) or (v3(A(x))[0]==0 and v3(A(x))[1]%3==2) or (v3(B(x))[0]%2==1 and B(x)%3**7!=0) or (v3(B(x))[0]==0 and v3(B(x))[1]%3==2))]
print("twisted 3-integral m surviving mod 3^7:",bad[:10],len(bad))
