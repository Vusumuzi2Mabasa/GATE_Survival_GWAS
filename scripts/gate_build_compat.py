"""Build-only fixes for upstream GATE 962f1028; preserve statistical routines."""
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('tools/GATE-main')
changes = {
    'NAMESPACE': [('useDynLib(SAIGE)', 'useDynLib(GATE)')],
    'R/RcppExports.R': [("PACKAGE = 'SAIGE'", "PACKAGE = 'GATE'")],
    'src/RcppExports.cpp': [('R_init_SAIGE(', 'R_init_GATE(')],
    'thirdParty/bgen/src/View.cpp': [('std::ios::streampos', 'std::streampos')],
    'configure': [('./waf\n', './waf -j1\n'),
                  ('rm -rf thirdParty/cget', 'mkdir -p thirdParty/cget'),
                  ('-DCMAKE_C_FLAGS="-fPIC"', '-DCMAKE_C_FLAGS="-fPIC -std=gnu11"'),
                  ('-fPIC -std=gnu11"', '-fPIC -std=gnu11 -Wno-error=implicit-function-declaration"')],
}
log = []
backup = Path('logs/gate/upstream_before_compat')
for p in (root / 'thirdParty/bgen').glob('.waf3-*/waflib/Node.py'):
    original = p.read_text()
    updated = original.replace('raise StopIteration', 'return')
    if updated != original:
        b = backup / 'waflib/Node.py'
        b.parent.mkdir(parents=True, exist_ok=True)
        if not b.exists():
            b.write_text(original)
        p.write_text(updated)
    log.append(str(p.relative_to(root)))
for name, replacements in changes.items():
    p = root / name
    original = p.read_text()
    updated = original
    for old, new in replacements:
        if old not in updated and new not in updated and not (
            name == 'configure' and '-fPIC -std=gnu11 -Wno-error=implicit-function-declaration' in updated
            and 'CMAKE_C_FLAGS' in old):
            raise RuntimeError(f'Unexpected source in {name}')
        updated = updated.replace(old, new)
    if updated != original:
        b = backup / name
        b.parent.mkdir(parents=True, exist_ok=True)
        if not b.exists():
            b.write_text(original)
        p.write_text(updated)
    log.append(name)
Path('logs/gate/build_compat.json').write_text(json.dumps({
    'source_commit': '962f1028120c73e24ec34f29e02eeae3bdc3bfb3',
    'files': log,
    'reason': 'Align GATE DLL lookup; use GNU C11 for legacy SuperLU; preserve completed dependencies; limit BGEN build jobs.',
    'statistical_algorithm_changes': False
}, indent=2))
