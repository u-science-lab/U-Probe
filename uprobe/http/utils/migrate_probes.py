"""One-time conversion: python -m uprobe.http.utils.migrate_probes SOURCE DESTINATION.

SOURCE is the original probe.json; DESTINATION is the new probes directory.
All original configuration fields survive exactly. Existing destinations are
never overwritten and the source file is left intact for deployment backup.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import tempfile

import yaml

from .probe_library import PresetLoader


def migrate(source, destination):
    original = json.loads(Path(source).read_text(encoding='utf-8-sig'))
    if not isinstance(original, dict) or not original:
        raise ValueError('Expected a non-empty method mapping')
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError(f'Refusing to overwrite {destination}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.probes-migration-', dir=destination.parent))
    try:
        for name, config in original.items():
            if not name or name in ('.', '..') or any(c in name for c in '/\\\0'):
                raise ValueError(f'Invalid method filename: {name!r}')
            if not isinstance(config, dict) or not isinstance(config.get('probes'), dict):
                raise ValueError(f'{name}: expected probes mapping')
            file = staging / f'{name}.yaml'
            file.write_text(yaml.safe_dump(config, sort_keys=False, allow_unicode=True), encoding='utf-8')
            if yaml.load(file.read_text(encoding='utf-8'), Loader=PresetLoader) != config:
                raise ValueError(f'{name}: YAML roundtrip changed the configuration')
        os.replace(staging, destination)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return list(original)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    print('Migrated methods:', ', '.join(migrate(args.source, args.destination)))
