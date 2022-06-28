import numpy as np
import sys
import pandas as pd
import itertools

n = int(sys.argv[1])  # how big of parameter sets we're checking
csvfile = sys.argv[2]  # where the ranking information is coming from
numrows = int(sys.argv[3])  # how many parameter sets are we using for combinations
numtosave = int(sys.argv[4])  # only going to save some subset of all the combinations checked
numcombos = int(np.prod([numrows - k for k in range(n)]) / np.prod(np.arange(n) + 1))  # ={numrows \choose n}
if numcombos < numtosave:
    numtosave = None
print(f'{numcombos} combinations will be checked.')  # so that if it's like a trillion then I just kill the script

stars, nums = ['HD1160', 'HR8799', 'KB', 'Kappa'], [13, 20, 30, 40, 50, 60]
columnnames = [f'{star}_{num}' for star in stars for num in nums]
df = pd.read_csv(csvfile)[:numrows]
indexcombos = itertools.combinations(np.arange(numrows), n)
data_to_save = {'Params': list(), 'Threshold': list()}
for ic in indexcombos:
    subdf = df.iloc[list(ic), :]
    threshold = np.min([np.max(subdf[col]) for col in columnnames])
    paramgroups = list()
    for idx in range(len(ic)):
        paramgroups.append(tuple([subdf.iloc[idx, k] for k in range(6)]))
    data_to_save['Params'].append(str(paramgroups))
    data_to_save['Threshold'].append(threshold)

broken_up_params = [list() for _ in range(n)]
for pset in data_to_save['Params']:
    pset = pset.replace("'", "")
    parens = list()
    for i in range(len(pset)):
        if pset[i] in ['(', ')']:
            parens.append(i)
    pts = [f'{pset[parens[2*j]: parens[2*j + 1] + 1]}' for j in range(int(len(parens) / 2))]
    assert len(pts) == len(broken_up_params)
    for k in range(n):
        broken_up_params[k].append(pts[k])
del data_to_save['Params']

combined = pd.DataFrame(data_to_save)
for m, psets in enumerate(broken_up_params):
    combined.insert(m, f'Params {m + 1}', psets)
to_csv = combined.sort_values('Threshold', ascending=False, ignore_index=True)[:numtosave]
to_csv.to_csv('groupscores_indsep.csv', index=False)
