import numpy as np
import sys
import pandas as pd
import itertools

n = int(sys.argv[1])  # how many parameter sets in each group we're checking
csvfile = sys.argv[2]  # where the ranking information is coming from (this should be output of parameterranking.py)
numrows = int(sys.argv[3])  # how many parameter sets are we using for combinations
numtosave = int(sys.argv[4])  # how many of the combinations checked to save in output file
output_filename = sys.argv[5]
numcombos = int(np.prod([numrows - k for k in range(n)]) / np.prod(np.arange(n) + 1))  # ={numrows \choose n}
if numcombos < numtosave or sys.argv[4].lower() == 'none':
    numtosave = None
print(f'{numcombos} combinations will be checked.')  # so that if it's like a trillion then I just kill the script

stars, nums = ['HD1160', 'HR8799', 'KB', 'Kappa'], [13, 20, 30, 40, 50, 60, 'SciSNR']
columnnames = [f'{star}_{num}' for star in stars for num in nums]
df = pd.read_csv(csvfile)[:numrows]
indexcombos = itertools.combinations(np.arange(numrows), n)


def find_threshold(ic):
    subdf = df.iloc[list(ic), :]
    threshold = np.min([np.max(subdf[col]) for col in columnnames])
    paramgroups = [tuple(subdf.iloc[idx, :5]) for idx in range(len(ic))]
    return str(paramgroups), threshold


data = list(map(find_threshold, indexcombos))
data_to_save = {'Params': [t[0] for t in data], 'Threshold': [t[1] for t in data]}


def break_up_params(pset):
    pset = pset.replace("'", "")
    parens = [i for i in range(len(pset)) if pset[i] in ['(', ')']]
    pts = [f'{pset[parens[2 * j]: parens[2 * j + 1] + 1]}' for j in range(int(len(parens) / 2))]
    return pts


bkp = list(map(break_up_params, data_to_save['Params']))
broken_up_params = [[pts[k] for pts in bkp] for k in range(n)]
del data_to_save['Params']

combined = pd.DataFrame(data_to_save)
for m, psets in enumerate(broken_up_params):
    combined.insert(m, f'Params {m + 1}', psets)
to_csv = combined.sort_values('Threshold', ascending=False, ignore_index=True)[:numtosave]
to_csv.to_csv(output_filename, index=False)
