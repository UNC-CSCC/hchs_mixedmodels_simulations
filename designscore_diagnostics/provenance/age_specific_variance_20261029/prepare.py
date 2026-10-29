"""Exact sufficient statistics for census Gaussian models; never accesses MI outputs."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import ctypes
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[3]  # Parent HCHS-Simulations repository, above this submodule.
if os.name == 'nt':
    kernel = ctypes.WinDLL('kernel32', use_last_error=True)
    kernel.GetCurrentProcess.restype = ctypes.c_void_p
    handle = kernel.GetCurrentProcess()
    kernel.SetPriorityClass.argtypes = (ctypes.c_void_p, ctypes.c_uint)
    kernel.SetProcessAffinityMask.argtypes = (ctypes.c_void_p, ctypes.c_size_t)
    assert kernel.SetPriorityClass(handle, 0x4000)
    assert kernel.SetProcessAffinityMask(handle, 1)

NAMES = ['(Intercept)', 'age_strat_new', 'x12', 'x17', 'x18', 'x6', 'y_bmi']

def r_vector(x):
    return 'c(' + ','.join(format(float(v), '.17g') for v in np.ravel(x)) + ')'

def r_matrix(x):
    return 'matrix(' + r_vector(x) + f',nrow={x.shape[0]},byrow=TRUE)'

def prepare(release):
    start = time.monotonic()
    source = ROOT / 'data' / release / f'population_3visits_{release}.csv'
    before = source.stat()
    cols = ['BGid', 'hhid', 'subid', 'v_num', 'age_base', 'x12', 'x17', 'x18', 'x6', 'y_bmi', 'y_gfr']
    print(f'{release}: reading complete population on one CPU at below-normal priority', flush=True)
    frame = pd.read_csv(source, usecols=cols)
    assert len(frame) % 3 == 0 and not frame.isna().any().any()
    n = len(frame) // 3
    assert np.array_equal(frame.v_num.to_numpy().reshape(n, 3), np.tile([1, 2, 3], (n, 1)))
    for name in ['BGid', 'hhid', 'subid', 'age_base']:
        values = frame[name].to_numpy().reshape(n, 3)
        assert np.all(values == values[:, :1]), name
    base = frame.iloc[::3]
    assert not base[['BGid', 'hhid', 'subid']].duplicated().any()
    household = pd.MultiIndex.from_frame(base[['BGid', 'hhid']])
    hcode, levels = pd.factorize(household, sort=False)
    assert np.all(hcode[1:] >= hcode[:-1])
    starts = np.r_[0, np.flatnonzero(hcode[1:] != hcode[:-1]) + 1]
    age = base.age_base.to_numpy() >= 45
    ng = np.bincount(age.astype(int), minlength=2)
    n0 = np.add.reduceat((~age).astype(int), starts)
    n1 = np.add.reduceat(age.astype(int), starts)
    counts, pattern = np.unique(np.column_stack([n0, n1]), axis=0, return_inverse=True)
    z = np.empty((n, 3, 8))
    z[:, :, 0] = 1
    z[:, :, 1] = age[:, None]
    for j, name in enumerate(NAMES[2:] + ['y_gfr'], start=2):
        z[:, :, j] = frame[name].to_numpy().reshape(n, 3)
    # Center every nonconstant column, including the outcome, for stable profiling.
    centers = z.mean(axis=(0, 1))
    centers[0] = 0
    z -= centers
    mean = z.mean(axis=1)
    subject_cross, within_cross, sums = [], [], []
    for g in range(2):
        take = age == bool(g)
        subject_cross.append(mean[take].T @ mean[take])
        contrast = (z[take] - mean[take, None, :]).reshape(-1, 8)
        within_cross.append(contrast.T @ contrast)
        sums.append(np.add.reduceat(mean * take[:, None], starts, axis=0))
    cross = []
    for i, (ny, no) in enumerate(counts):
        take = pattern == i
        u0, u1 = sums[0][take], sums[1][take]
        cross.append(dict(n0=int(ny), n1=int(no), households=int(take.sum()),
                          c00=u0.T @ u0, c11=u1.T @ u1, c01=u0.T @ u1))
    fields = [f'n={n}', f'nh={len(starts)}', 'ng=' + r_vector(ng), 'center=' + r_vector(centers),
              'mean=list(' + ','.join(map(r_matrix, subject_cross)) + ')',
              'within=list(' + ','.join(map(r_matrix, within_cross)) + ')']
    blocks = []
    for item in cross:
        blocks.append('list(' + ','.join(f'{k}={r_matrix(v) if isinstance(v, np.ndarray) else v}'
                                        for k, v in item.items()) + ')')
    fields.append('patterns=list(' + ','.join(blocks) + ')')
    (OUT / f'sufficient_{release}.R').write_text('list(' + ',\n'.join(fields) + ')\n')
    # Small observed-data subset for an independent direct covariance likelihood check.
    test = frame.iloc[: min(120, len(frame))].copy()
    # R's read.table can parse this small tabular validation input using base R.
    test.to_csv(OUT / f'validation_subset_{release}.txt', sep='\t', index=False)
    after = source.stat()
    assert before.st_size == after.st_size and before.st_mtime_ns == after.st_mtime_ns
    record = dict(population=release, subjects=n, households=len(starts),
                  younger=int(ng[0]), older=int(ng[1]), household_patterns=len(counts),
                  source=str(source), source_bytes=before.st_size, source_mtime_ns=before.st_mtime_ns,
                  cpu_affinity='one logical CPU', process_priority='below normal',
                  fixed_effects=NAMES, elapsed_seconds=time.monotonic()-start)
    (OUT / f'input_{release}.json').write_text(json.dumps(record, indent=2))
    print(json.dumps(record), flush=True)

if __name__ == '__main__':
    for release in ['Dec2024', 'Oct2026']:
        prepare(release)
