"""Build the household lookup for scen3_designscore without fitting any model.

Reads identifiers, visit, age, stratum and surname-group membership only.
The input must be the full saved population, ordered by subject and visit.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

BG_COUNTS = (58, 21, 130, 167) * 2
BG_PROBS = tuple(x / 1.8 for x in (.25, .25, .60, .60)) * 2
HH_HISP = (.18, .225, .14, .14) * 2
HH_OTHER = (.025, .0175, .035, .04) * 2
FIELDS = ('population', 'bgid', 'hhid', 'strat', 'm_h', 'p_h',
          'n_younger', 'n_older')


def integer(value: str) -> int:
    number = float(value)
    if not math.isfinite(number) or number != math.floor(number):
        raise ValueError(f'Expected an integer, received {value!r}')
    return int(number)


def logical(value: str) -> bool:
    value = value.strip().lower()
    if value in ('true', '1'):
        return True
    if value in ('false', '0'):
        return False
    raise ValueError(f'Invalid hisp_strat value {value!r}')


def read_households(path: Path) -> dict:
    households = {}
    previous_key = None
    expected_visit = 1
    person_info = None
    with path.open(newline='', encoding='utf-8-sig') as source:
        reader = csv.DictReader(source)
        names = {name.lower(): name for name in (reader.fieldnames or [])}
        required = ('bgid', 'hhid', 'subid', 'strat', 'v_num', 'hisp_strat', 'age_base')
        if not set(required).issubset(names):
            raise ValueError(f'Missing population columns: {set(required) - names.keys()}')
        for row in reader:
            bg, hh, sub, strat, visit = (integer(row[names[k]]) for k in
                                        ('bgid', 'hhid', 'subid', 'strat', 'v_num'))
            age = float(row[names['age_base']])
            hisp = logical(row[names['hisp_strat']])
            if min(bg, hh, sub) < 1 or not 1 <= strat <= 8 or not math.isfinite(age):
                raise ValueError('Invalid population design identifiers or age')
            key = (bg, hh, sub)
            info = (strat, hisp, age)
            if expected_visit == 1:
                if visit != 1 or (previous_key is not None and key <= previous_key):
                    raise ValueError('Require unique subjects sorted by BG, household, subject, visit')
                person_info = info
                household = households.setdefault((bg, hh), {
                    'strat': strat, 'hisp': hisp, 'm_h': 0, 'n_younger': 0, 'n_older': 0})
                if (household['strat'], household['hisp']) != (strat, hisp):
                    raise ValueError('Household spans strata or surname groups')
                household['m_h'] += 1
                household['n_older' if age >= 45 else 'n_younger'] += 1
            elif key != previous_key or info != person_info or visit != expected_visit:
                raise ValueError('Incomplete or inconsistent three-visit subject')
            previous_key = key
            expected_visit = 1 if visit == 3 else visit + 1
    if not households or expected_visit != 1:
        raise ValueError('Empty population or incomplete last subject')
    return households


def frame_rows(households: dict, population: str) -> list[dict]:
    bg_stratum = {}
    counts = Counter()
    for (bg, hh), row in households.items():
        if bg_stratum.setdefault(bg, row['strat']) != row['strat']:
            raise ValueError('Block group spans strata')
        counts[bg, row['hisp']] += 1
    observed = Counter(bg_stratum.values())
    if tuple(observed[s] for s in range(1, 9)) != BG_COUNTS:
        raise ValueError('Population BG counts do not match the HCHS sampling configuration')
    result = []
    for (bg, hh), row in sorted(households.items()):
        s = row['strat'] - 1
        n_group = counts[bg, row['hisp']]
        fraction = HH_HISP[s] if row['hisp'] else HH_OTHER[s]
        # Python and the production R generator both round ties to even.
        p_bg = round(BG_PROBS[s] * BG_COUNTS[s]) / BG_COUNTS[s]
        p_hh = round(fraction * n_group) / n_group
        if not 0 < p_bg * p_hh <= 1:
            raise ValueError(f'Zero/invalid household inclusion probability for {(bg, hh)}')
        result.append(dict(population=population, bgid=bg, hhid=hh, strat=s + 1,
                           m_h=row['m_h'], p_h=p_bg * p_hh,
                           n_younger=row['n_younger'], n_older=row['n_older']))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--population', choices=('dec2024', 'oct2026'), required=True)
    parser.add_argument('--population-file', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    manifest = output.with_suffix('.json')
    if output.exists() or manifest.exists():
        parser.error('Output or provenance file already exists; choose a new path')
    households = read_households(args.population_file)
    rows = frame_rows(households, args.population)
    with args.population_file.open('rb') as source:
        digest = hashlib.file_digest(source, 'sha256').hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', newline='', encoding='utf-8') as target:
        writer = csv.DictWriter(target, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    with manifest.open('x', encoding='utf-8') as target:
        json.dump(dict(population=args.population, source=str(args.population_file.resolve()),
                       source_sha256=digest, households=len(rows),
                       subjects=sum(r['m_h'] for r in rows),
                       design='HCHS stage 1/2 rounded SRS fractions; no outcome fields used'),
                  target, indent=2)
    print(f'Wrote {len(rows):,} households to {output}. No models were fitted.')


if __name__ == '__main__':
    main()
