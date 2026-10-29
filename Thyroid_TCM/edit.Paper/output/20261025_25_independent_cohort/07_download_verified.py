from pathlib import Path
import requests,concurrent.futures,gzip,hashlib,json,time,urllib3
urllib3.disable_warnings();p=Path(__file__).parent
u='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33630/matrix/GSE33630_series_matrix.txt.gz'
N=31398276;step=1048576
parts=p/'08_download_parts';parts.mkdir(exist_ok=True)
def get(start):
 end=min(start+step,N)-1;f=parts/str(start)
 if f.exists() and f.stat().st_size==end-start+1:return
 for attempt in range(4):
  try:
   r=requests.get(u+f'?part={start}',headers={'Range':f'bytes={start}-{end}'},timeout=(15,45),verify=False)
   assert r.status_code==206,(r.status_code,start)
   assert r.headers.get('Content-Range')==f'bytes {start}-{end}/{N}',r.headers.get('Content-Range')
   assert len(r.content)==end-start+1
   f.write_bytes(r.content);print('part',start,'ok',flush=True);return
  except Exception as e: print('retry',start,attempt,str(e)[:100],flush=True)
 raise RuntimeError(start)
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:list(ex.map(get,range(0,N,step)))
b=b''.join((parts/str(i)).read_bytes() for i in range(0,N,step));assert len(b)==N;gzip.decompress(b)
f=p/'02_GSE33630_series_matrix.txt.gz';f.write_bytes(b)
(p/'03_download_manifest.json').write_text(json.dumps({'url':u,'bytes':N,'sha256':hashlib.sha256(b).hexdigest()},indent=2));print('VERIFIED COMPLETE',flush=True)
r=requests.get('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33630/suppl/GSE33630_clinical-annotation.txt.gz',timeout=(15,40),verify=False);r.raise_for_status();(p/'06_clinical_annotation.txt').write_bytes(gzip.decompress(r.content));print('clinical annotation saved',flush=True)
