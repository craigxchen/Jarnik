from f27 import *
labs=[(x,y) for x in range(27) for y in range(27)]
FORMS={'x':lambda x,y:x,'y':lambda x,y:y,'x+y':lambda x,y:ADD[x][y],'x-y':lambda x,y:ADD[x][neg(y)]}
def rows(fs,sg):
    return [[1+CHI[FORMS[f](x,y) if s==1 else neg(FORMS[f](x,y))] for (x,y) in labs] for f,s in zip(fs,sg)]
names=list(FORMS)
for m in (2,3,4):
    for fs in itertools.combinations(names,m):
        res=[]
        for sg in itertools.product((1,-1),repeat=m):
            if sg[0]==-1: continue
            R=rows(fs,sg); g=sum(min(c) for c in zip(*R))
            res.append((g,sg))
        print(m,fs,sorted(res,reverse=True)[:3])
