"""Download PubMed XML for the archived search; never overwrite source files."""
import json, pathlib, subprocess, datetime, hashlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / 'output/01_检索与研究方案/raw'
search = json.loads((RAW / 'pubmed_esearch.json').read_text())['esearchresult']
ids = search['idlist']
assert len(ids) == int(search['count']) == len(set(ids))
for start in range(0, len(ids), 150):
    dest = RAW / f'pubmed_{start:04d}.xml'
    if dest.exists():
        print('Already downloaded:', dest.name)
        continue
    subprocess.run(['curl', '-sS', '-L', '--fail', '--retry', '3', '--max-time', '90', '-G',
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi', '--data-urlencode', 'db=pubmed',
        '--data-urlencode', 'retmode=xml', '--data-urlencode', 'id=' + ','.join(ids[start:start+150]),
        '-o', str(dest)], check=True)
    print('Downloaded:', dest.name, dest.stat().st_size, flush=True)
manifest = {'retrieved_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'database':'PubMed', 'count':len(ids), 'query_translation':search['querytranslation'],
    'publication_cutoff':'2026-10-10', 'language_filter_at_search':False,
    'note':'Publication-date filter is not a historical database-indexing snapshot.',
    'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in RAW.glob('*') if p.name != 'manifest.json'}}
(RAW / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
