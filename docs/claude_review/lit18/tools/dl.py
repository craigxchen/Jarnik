import sys, os, urllib.request, concurrent.futures, ssl
cat = sys.argv[1]; outdir = sys.argv[2]; dates = open(sys.argv[3]).read().split()
sha = '30c967b1cf20a5056694b79f6ad1bca8b3e7cf89'
os.makedirs(outdir, exist_ok=True)
ctx = ssl.create_default_context(cafile=os.environ.get('SSL_CERT_FILE') or '/root/.ccr/ca-bundle.crt')
def get(d):
    p = os.path.join(outdir, d + '.xml')
    if os.path.exists(p): return 'have'
    url = f'https://raw.githubusercontent.com/ehijano/rss_fetch/{sha}/rss_data/{cat}/{d}_{cat}.xml'
    try:
        with urllib.request.urlopen(url, timeout=40, context=ctx) as r:
            data = r.read()
        open(p, 'wb').write(data); return 'ok'
    except Exception as e:
        return 'miss'
with concurrent.futures.ThreadPoolExecutor(12) as ex:
    res = list(ex.map(get, dates))
print(cat, {k: res.count(k) for k in set(res)})
