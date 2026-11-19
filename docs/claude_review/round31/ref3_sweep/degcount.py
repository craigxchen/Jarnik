# arithmetic of the fifteen-cut and ten-pair degree cutoffs (III.10)
import math
for d in range(1,7):
    D=15*d; k=math.ceil(D/2); contact=8*d-d//2
    print(d,'D',D,'k=ceil(D/2)',k,'8d-floor(d/2)',contact,'ok',contact>=k,'4k+1>2D',4*k+1>2*D,'4k<=30d+2',4*k<=30*d+2)
# ten-pair: projected 3x3 det poles 18d, order 3*6d+1
for d in range(1,5): print('tenpair d',d,'poles',18*d,'order',18*d+1)
print('ten-pair cutoff c<1/(24 log 6)=',1/(24*math.log(6)),' fifteen-cut: degree m^{c*30*log 8}, o(m) iff c<',1/(30*math.log(8)))
