from pathlib import Path
import requests,gzip,hashlib,json
p=Path(__file__).parent
s=requests.Session();s.trust_env=False
u='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33630/matrix/GSE33630_series_matrix.txt.gz'
r=s.get(u,timeout=(20,60));r.raise_for_status();gzip.decompress(r.content)
f=p/'02_GSE33630_series_matrix.txt.gz';f.write_bytes(r.content)
(p/'03_download_manifest.json').write_text(json.dumps({'url':u,'bytes':len(r.content),'sha256':hashlib.sha256(r.content).hexdigest()},indent=2))
t=gzip.decompress(r.content).decode()
for line in t.splitlines():
 if line.startswith(('!Series_title','!Series_summary','!Series_overall_design','!Series_pubmed_id','!Sample_title','!Sample_characteristics_ch1','!Sample_data_processing')): print(line[:1800])
