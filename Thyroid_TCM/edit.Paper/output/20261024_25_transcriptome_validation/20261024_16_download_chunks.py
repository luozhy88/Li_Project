import requests,concurrent.futures,time,gzip,os
from pathlib import Path
base=Path('edit.Paper/output/20261024_25_transcriptome_validation')
items=[('20261024_GSE60542_series_matrix.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60542/matrix/GSE60542_series_matrix.txt.gz'),('20261024_GPL570.annot.gz','https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz')]
# Discover exact object lengths via a one-byte range request; never rely on estimates.
def remote_size(url):
 s=requests.Session();s.trust_env=False
 r=s.get(url,headers={'Range':'bytes=0-0'},timeout=(20,60),stream=True)
 if r.status_code!=206 or 'Content-Range' not in r.headers: raise RuntimeError(f'Cannot verify remote size: {r.status_code} {r.headers}')
 total=int(r.headers['Content-Range'].split('/')[-1]);r.close();return total
items=[(name,url,remote_size(url)) for name,url in items]
print('VERIFIED_REMOTE_LENGTHS',items,flush=True)
# Each chunk verifies the returned byte range before being combined.
def fetch(job):
 name,url,total,start,end=job;p=base/(name+f'.part{start}')
 if p.exists() and p.stat().st_size==end-start+1:return
 for k in range(5):
  try:
   s=requests.Session();s.trust_env=False
   r=s.get(url,headers={'Range':f'bytes={start}-{end}'},timeout=(20,60))
   assert r.status_code==206,(r.status_code,r.headers)
   cr=r.headers['Content-Range']; actual_total=int(cr.split('/')[-1])
   assert cr.startswith(f'bytes {start}-{min(end,actual_total-1)}/'),cr
   assert len(r.content)==min(end,actual_total-1)-start+1
   p.write_bytes(r.content);print(name,start,len(r.content),actual_total,flush=True);return
  except Exception as e:
   print('retry',name,start,k,str(e)[:120],flush=True);time.sleep(1)
 raise RuntimeError(str(job))
jobs=[]
for name,url,total in items:
 for a in range(0,total,1024*1024): jobs.append((name,url,total,a,min(a+1024*1024-1,total-1)))
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:list(ex.map(fetch,jobs))
for name,url,total in items:
 tmp=base/(name+'.complete')
 with tmp.open('wb') as out:
  for a in range(0,total,1024*1024):out.write((base/(name+f'.part{a}')).read_bytes())
 with gzip.open(tmp,'rb') as f:
  while f.read(1024*1024):pass
 print('COMPLETE',name,tmp.stat().st_size,flush=True)
