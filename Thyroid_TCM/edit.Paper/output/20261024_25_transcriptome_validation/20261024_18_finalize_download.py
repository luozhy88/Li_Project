import gzip
import hashlib
import os
from pathlib import Path

base = Path(__file__).parent
items = [
    ("20261024_GSE60542_series_matrix.txt.gz", 27530495),
    ("20261024_GPL570.annot.gz", 8471521),
]

for name, total in items:
    complete = base / f"{name}.complete"
    with complete.open("wb") as out:
        for start in range(0, total, 1048576):
            end = min(start + 1048576, total) - 1
            part = base / f"{name}.part{start}"
            if not part.is_file() or part.stat().st_size != end - start + 1:
                raise SystemExit(f"missing or short {part}")
            with part.open("rb") as inp:
                while chunk := inp.read(1048576):
                    out.write(chunk)
    if complete.stat().st_size != total:
        raise SystemExit(f"wrong size {complete}")
    with gzip.open(complete, "rb") as inp:
        while inp.read(1048576):
            pass
    digest = hashlib.sha256(complete.read_bytes()).hexdigest()
    os.replace(complete, base / name)
    print("VALID", name, total, digest, flush=True)
