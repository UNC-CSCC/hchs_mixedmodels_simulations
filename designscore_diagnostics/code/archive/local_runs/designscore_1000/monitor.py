import csv
import json
from collections import Counter
from datetime import datetime
from pathlib import Path

state = Path('.statatmp/designscore_1000')
progress = json.loads((state/'status.json').read_text())
folder = Path('data/estimates/oct2026/scen3_designscore/exact')
fits = list(folder.glob('scen3_designscore_coef_*.csv'))
status = Counter()
for path in fits:
    sim = path.stem.rsplit('_', 1)[1]
    with (folder/f'scen3_designscore_audit_{sim}.csv').open() as source:
        audit = next(csv.DictReader(source))
        status[audit['fit_status']] += 1
stages = Counter(job['stage'] for job in progress['jobs'].values())
errors = {label: job for label,job in progress['jobs'].items()
          if job['stage']=='error' or job.get('failed_fits')}
manifest_failures = []
for first in range(201, 1001, 50):
    path = folder/f'local_run_{first}_{first+49}.csv'
    if path.exists():
        with path.open() as source:
            manifest_failures.extend(row for row in csv.DictReader(source)
                if row.get('status') == 'failed')
print(json.dumps(dict(time=datetime.now().isoformat(timespec='seconds'),
    status=progress['status'], completed=len(fits), new_completed=len(fits)-200,
    fit_status=dict(status), job_stages=dict(stages), errors=errors,
    manifest_failures=manifest_failures)))
