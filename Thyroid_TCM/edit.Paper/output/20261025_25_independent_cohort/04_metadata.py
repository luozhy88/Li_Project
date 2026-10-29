from pathlib import Path
import requests,zlib,gzip
p=Path(__file__).parent;s=requests.Session();s.trust_env=False
u='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33630/'
r=s.get(u+'matrix/GSE33630_series_matrix.txt.gz',headers={'Range':'bytes=0-131071'},timeout=(15,40));r.raise_for_status();print(r.status_code,r.headers.get('Content-Range'),flush=True)
t=zlib.decompressobj(16+zlib.MAX_WBITS).decompress(r.content).decode()
meta='\n'.join(x for x in t.splitlines() if x.startswith('!'))
(p/'05_series_metadata.txt').write_text(meta)
for line in meta.splitlines():
 if line.startswith(('!Series_title','!Series_summary','!Series_overall_design','!Series_pubmed_id','!Sample_title','!Sample_characteristics_ch1','!Sample_data_processing')):print(line[:2500],flush=True)
r=s.get(u+'suppl/GSE33630_clinical-annotation.txt.gz',timeout=(15,40));r.raise_for_status();(p/'06_clinical_annotation.txt').write_bytes(gzip.decompress(r.content));print(gzip.decompress(r.content).decode()[:2500])
