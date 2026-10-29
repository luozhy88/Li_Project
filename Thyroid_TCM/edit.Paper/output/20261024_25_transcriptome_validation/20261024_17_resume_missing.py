import os
import time
from pathlib import Path

import requests

base = Path(__file__).parent
jobs = [
    (
        "20261024_GSE60542_series_matrix.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60542/matrix/GSE60542_series_matrix.txt.gz",
        14680064,
        15728639,
        27530495,
    ),
    (
        "20261024_GPL570.annot.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
        7340032,
        8388607,
        8471521,
    ),
]

for name, url, start, end, total in jobs:
    target = base / f"{name}.part{start}"
    temp = target.with_suffix(target.suffix + ".resume")
    ok = False
    for attempt in range(30):
        try:
            data = bytearray()
            pos = start
            while pos <= end:
                sub_end = min(pos + 262143, end)
                r = requests.get(
                    url,
                    headers={"Range": f"bytes={pos}-{sub_end}"},
                    timeout=(30, 90),
                    proxies={"http": None, "https": None},
                )
                if r.status_code != 206:
                    raise RuntimeError(f"HTTP {r.status_code}")
                cr = r.headers.get("Content-Range", "")
                expected = sub_end - pos + 1
                if not cr.startswith(f"bytes {pos}-{sub_end}/") or len(r.content) != expected:
                    raise RuntimeError(f"short range {cr} len={len(r.content)} expected={expected}")
                data.extend(r.content)
                pos = sub_end + 1
            if len(data) != end - start + 1:
                raise RuntimeError(f"wrong total {len(data)}")
            temp.write_bytes(data)
            os.replace(temp, target)
            print("DONE", name, start, len(data), total, flush=True)
            ok = True
            break
        except Exception as exc:
            print("RETRY", name, attempt + 1, str(exc)[:160], flush=True)
            time.sleep(2)
    if not ok:
        raise SystemExit(f"FAILED {name} {start}-{end}")
