import sys, os, re, json, html
root = sys.argv[1]; out = sys.argv[2]
db = {}
dates = []
for cat in sorted(os.listdir(root)):
    p = os.path.join(root, cat)
    if not os.path.isdir(p): continue
    for f in sorted(os.listdir(p)):
        dates.append(f[:10])
        s = open(os.path.join(p, f), encoding='utf-8', errors='replace').read()
        for m in re.finditer(r'<item>(.*?)</item>', s, re.S):
            it = m.group(1)
            t = re.search(r'<title>(.*?)</title>', it, re.S)
            l = re.search(r'<link>(.*?)</link>', it, re.S)
            d = re.search(r'<description>(.*?)</description>', it, re.S)
            a = re.search(r'<dc:creator>(.*?)</dc:creator>', it, re.S)
            if not (t and l): continue
            aid = l.group(1).strip().split('/')[-1]
            desc = html.unescape(d.group(1)) if d else ''
            ver = re.search(r'arXiv:\S+?v(\d+)', desc)
            ver = int(ver.group(1)) if ver else 0
            rec = dict(id=aid, title=html.unescape(t.group(1)).strip(), authors=html.unescape(a.group(1)).strip() if a else '', abstract=re.sub(r'^arXiv:\S+\s+Announce Type:\s*\S+\s*Abstract:\s*', '', desc.strip()), date=f[:10], cat=cat, ver=ver)
            if aid not in db or ver >= db[aid]['ver']:
                first = db.get(aid, {}).get('first', f[:10])
                rec['first'] = min(first, f[:10]); db[aid] = rec
json.dump(db, open(out, 'w'))
print(len(db), min(dates), max(dates))
