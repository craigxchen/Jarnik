import json,re,sys
db=json.load(open(sys.argv[1]))
pats=sys.argv[2:]
hits=[]
for k,r in db.items():
    t=(r['title']+' '+r['abstract'])
    if all(re.search(p,t,re.I|re.S) for p in pats):
        hits.append(r)
hits.sort(key=lambda r:r['id'])
print(len(hits),'hits for',pats)
for r in hits[:60]:
    print(r['id'],'|',r['title'][:110],'|',r['authors'][:60],'|',r['cat'])
