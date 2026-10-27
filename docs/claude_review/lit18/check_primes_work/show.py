import sys, json
db = json.load(open(sys.argv[1]))
for i in sys.argv[2:]:
    r = db.get(i)
    if not r: print('MISSING-FROM-MIRROR', i, '\n'); continue
    print(f"### {r['id']} [first {r['first']} .. last {r['date']} v{r['ver']} {r['cat']}] {r['title']}\n   AUTH: {r['authors'][:300]}\n   ABS: {r['abstract'][:3000]}\n")
