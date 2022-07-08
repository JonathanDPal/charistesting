import numpy as np
import sys
import pandas as pd
import itertools
from copy import copy

N = int(sys.argv[1])  # how big of parameter sets we're checking
csvfile = sys.argv[2]  # where the ranking information is coming from
output_filename = sys.argv[3]
num_to_skip = int(sys.argv[4])  # ones that have already been checked
max_threshold = float(sys.argv[5]) if '.csv' not in sys.argv[5] else sys.argv[5]  # start at current known value
if len(sys.argv) == 8:
    batched = True
    batchindex, batchsize = int(sys.argv[6]), int(sys.argv[7])
else:
    batched = False

if isinstance(max_threshold, str):
    max_threshold = pd.read_csv(max_threshold)['Threshold'][0]

stars, nums = ['HD1160', 'HR8799', 'KB', 'Kappa'], [13, 20, 30, 40, 50, 60, 'SciSNR']
columnnames = [f'{star}_{num}' for star in stars for num in nums]
df = pd.read_csv(csvfile)
numrows = len(df)
indexcombos = [ic for ic in itertools.combinations(np.arange(numrows), N-1)][num_to_skip: None]
numindexcombos = len(indexcombos)
if batched:
    remainder = numindexcombos % batchsize
    partial_batch = False
    if remainder != 0:
        num_full_size_batches = int(np.floor(numindexcombos / batchsize))
        if batchindex > num_full_size_batches:
            partial_batch = True
    startingindex = batchsize * (batchindex - 1)  # (batchindex - 1) is because of 1-based indexing
    if not partial_batch:
        finalindex = startingindex + batchsize
    else:
        finalindex = None
    indexcombos = indexcombos[startingindex: finalindex]
    numindexcombos = len(indexcombos)

if not batched:  # don't want a bunch of print statements on the batched stuff
    print(f'{numindexcombos} combinations will be checked.')
data_to_save = {'Params': list(), 'Threshold': list()}

for idx, ic in enumerate(indexcombos):
    ic = list(ic)
    subdf = df.iloc[ic, :]
    if subdf['max'].min() < max_threshold:
        continue
    colmaxes = pd.DataFrame({'max': [np.max(subdf[col]) for col in columnnames]}, index=columnnames).sort_values('max')
    curr = np.max(ic) + 1
    for k in np.arange(numrows - curr) + curr:
        found = 0
        if found < N - 1 and k in ic:
            found += 1
            continue
        ssdf = df.iloc[k, :]
        if ssdf['max'] < max_threshold:
            continue
        local_cm = copy(colmaxes)
        threshold = local_cm['max'][0]
        for col, row in local_cm.iterrows():
            new_val = ssdf[col]
            old_val = row['max']
            if new_val < old_val:
                if old_val < max_threshold:
                    break
            else:
                if new_val < max_threshold:
                    break
                else:
                    row['max'] = new_val
                    threshold = local_cm['max'].min()
        if threshold > max_threshold:
            ssubdf = copy(subdf)
            ssubdf.loc[len(ssubdf.index)] = ssdf
            ssubdf.index = np.arange(N)
            paramgroups = [tuple(ssubdf.iloc[idx, :5]) for idx in range(N)]
            data_to_save['Params'] = [str(paramgroups)]
            data_to_save['Threshold'] = [threshold]
            max_threshold = threshold
        elif threshold == max_threshold:
            ssubdf = copy(subdf)
            ssubdf.loc[len(ssubdf.index)] = ssdf
            ssubdf.index = np.arange(N)
            paramgroups = [tuple(ssubdf.iloc[idx, :5]) for idx in range(N)]
            data_to_save['Params'].append(str(paramgroups))
            data_to_save['Threshold'].append(threshold)


def break_up_params(pset):
    pset = pset.replace("'", "")
    parens = [i for i in range(len(pset)) if pset[i] in ['(', ')']]
    pts = [f'{pset[parens[2 * j]: parens[2 * j + 1] + 1]}' for j in range(int(len(parens) / 2))]
    return pts


bkp = list(map(break_up_params, data_to_save['Params']))
broken_up_params = [[pts[k] for pts in bkp] for k in range(N)]
del data_to_save['Params']
combined = pd.DataFrame(data_to_save)
for m, psets in enumerate(broken_up_params):
    combined.insert(m, f'Params {m + 1}', psets)
to_save = combined.sort_values('Threshold', ascending=False, ignore_index=True)
to_save.to_csv(output_filename, index=False)
