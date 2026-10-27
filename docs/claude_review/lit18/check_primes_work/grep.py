import sys, json, re
db = json.load(open(sys.argv[1]))
maxn = int(sys.argv[2])
pats = [re.compile(p, re.I|re.S) for p in sys.argv[3:]]
n=0
for k, r in sorted(db.items(), key=lambda kv: kv[1]['first']):
    txt = r['title'] + ' ' + r['abstract']
    if all(p.search(txt) for p in pats):
        n+=1
        if n<=maxn:
            print(f"### {r['id']} [{r['first']}..{r['date']} v{r['ver']} {r['cat']}] {r['title']}\n   AUTH: {r['authors'][:120]}\n   {r['abstract'][:500]}\n")
print('matches', n)
