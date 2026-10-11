from pathlib import Path
import zipfile
import tarfile
import json

bad = []
root = Path('/mnt/gate-runtime/cache')
for p in root.iterdir():
    try:
        if p.name.endswith('.conda'):
            with zipfile.ZipFile(p) as z:
                assert z.testzip() is None
        elif p.name.endswith('.tar.bz2'):
            with tarfile.open(p) as t:
                t.getmembers()
    except Exception as e:
        bad.append({'file': p.name, 'bytes': p.stat().st_size, 'error': str(e)})
print(json.dumps(bad, indent=2), flush=True)
Path('logs/gate/cache-integrity.json').write_text(json.dumps(bad, indent=2))
