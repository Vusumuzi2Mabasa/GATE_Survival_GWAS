"""Install the cached, exact environment inside the Toshiba-backed Linux image."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import os
import zipfile

project = Path('/mnt/d/survival_gwas_simulation')
runtime = Path('/mnt/gate-runtime')
subprocess.run(['mountpoint', '-q', str(runtime)], check=True)
cache = runtime / 'cache'
cache.mkdir(exist_ok=True)
records = []
archives = (project / 'logs/gate/cache-archives.txt').read_text().splitlines()
for k, name in enumerate(archives):
    src = project / name
    dst = cache / src.name
    if not dst.exists() or dst.stat().st_size != src.stat().st_size:
        shutil.copyfile(src, dst)
    content = dst.read_bytes()
    records.append({'file': dst.name, 'sha256': hashlib.sha256(content).hexdigest(),
                    'md5': hashlib.md5(content).hexdigest()})
    if k % 20 == 0:
        print('Cached archives copied:', k, flush=True)
explicit = runtime / 'environment-explicit-local.txt'
explicit.write_text('@EXPLICIT\n' + '\n'.join((cache / r['file']).as_uri() + '#' + r['md5'] for r in records) + '\n')
(project / 'logs/gate/archive-checksums.json').write_text(json.dumps(records, indent=2))
env = os.environ.copy()
env['MAMBA_ROOT_PREFIX'] = str(runtime / 'mamba')
env['MAMBA_EXTRACT_THREADS'] = '2'
mm = project / 'tools/micromamba/bin/micromamba'
subprocess.run([str(mm), 'create', '-y', '-p', str(runtime / 'env'), '--file', str(explicit)], env=env, check=True)
# Source is extracted from the preserved upstream archive, not the edited copy.
source = runtime / 'source'
source.mkdir(exist_ok=True)
with zipfile.ZipFile(project / 'tools/GATE-main.zip') as z:
    for i in z.infolist():
        n = i.filename
        if i.is_dir() or '/.hg/' in n or '/thirdParty/cget/' in n or '/thirdParty/bgen/build/' in n:
            continue
        if '/boost_1_55_0/' in n:
            b = n.split('/boost_1_55_0/')[1]
            if not (b.startswith('boost/') or b == 'wscript' or (b.startswith('libs/') and '/src/' in b)):
                continue
        z.extract(i, source)
print('Isolated environment and source ready inside Toshiba image.', flush=True)
