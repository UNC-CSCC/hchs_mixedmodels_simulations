import csv
import json
from collections import Counter
from datetime import datetime
from pathlib import Path

root = Path('data/estimates/oct2026/scen3_designscore/exact')
fits = list(root.glob('scen3_designscore_coef_*.csv'))
status = Counter()
for path in fits:
    sim = path.stem.rsplit('_', 1)[1]
    with (root/f'scen3_designscore_audit_{sim}.csv').open() as source:
        row = next(csv.DictReader(source))
        status[row['fit_status']] += 1
ranges, errors = {}, []
for first in range(1, 201, 25):
    path = root/f'local_run_{first}_{first+24}.csv'
    rows = list(csv.DictReader(path.open())) if path.exists() else []
    ranges[f'{first}-{first+24}'] = len(rows)
    errors.extend(row for row in rows if row.get('return_code') not in ('0', None))
print(json.dumps(dict(time=datetime.now().isoformat(timespec='seconds'),
    completed=len(fits), fit_status=dict(status), ranges=ranges, errors=errors)))
