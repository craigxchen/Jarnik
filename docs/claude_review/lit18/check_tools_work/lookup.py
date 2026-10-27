import sys, json, re
db = json.load(open(sys.argv[1]))
ids = [a for a in sys.argv[2:] if not a.startswith('re:')]
pats = [a[3:] for a in sys.argv[2:] if a.startswith('re:')]
for i in ids:
    r = db.get(i)
    if not r: print('MISSING', i); continue
    print(f"### {r['id']} [{r['first']}..{r['date']} v{r['ver']} {r['cat']}] {r['title']}\n   AUTH: {r['authors'][:200]}\n   {r['abstract'][:1500]}\n")
if pats:
    ps=[re.compile(p, re.I) for p in pats]; n=0
    for k,r in sorted(db.items(), key=lambda kv: kv[1]['first']):
        t=r['title']+' '+r['abstract']+' '+r['authors']
        if all(p.search(t) for p in ps):
            n+=1; print(f"* {r['id']} [{r['first']}..{r['date']} v{r['ver']} {r['cat']}] {r['title']} | {r['authors'][:120]}")
    print('matches',n)
