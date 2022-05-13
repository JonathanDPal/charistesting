import numpy as np
import sys
import pandas as pd
import itertools

n = int(sys.argv[1])  # how big of parameter sets we're checking
csvfile = sys.argv[2]  # where the ranking information is coming from
numrows = int(sys.argv[3])  # how many parameter sets are we using for combinations
numtosave = int(sys.argv[4])  # only going to save some subset of all the combinations checked
if len(sys.argv) == 6:  # can specify to only look at contrast/snr scores
    if str.lower(sys.argv[5]) == 'contrast':
        contrastonly, snronly = True, False
    elif str.lower(sys.argv[5]) == 'snr':
        contrastonly, snronly = False, True
elif len(sys.argv) == 7:  # can specify to only do Kappa/exclude Kappa
    if str.lower(sys.argv[5]) == 'contrast':
        contrastonly, snronly = True, False
    elif str.lower(sys.argv[5]) == 'snr':
        contrastonly, snronly = False, True
    else:
        contrastonly, snronly = False, False
    if str.lower(sys.argv[6]) in ['b', 'broadband']:
        broadbandonly, kappaonly = True, False
    elif str.lower(sys.argv[6]) in ['k', 'kappa']:
        broadbandonly, kappaonly = False, True
    else:
        broadbandonly, kappaonly = False, False
else:
    contrastonly, snronly, broadbandonly, kappaonly = False, False, False, False
numcombos = int(np.prod([numrows - k for k in range(n)]) / np.prod(np.arange(n) + 1))  # ={numrows \choose n}
print(f'{numcombos} combinations will be checked.')  # so that if it's like a trillion then I just kill the script

if contrastonly:
    if broadbandonly:
        columnnames = ['HD Contrast', 'HR Contrast', 'HIP Contrast']
    elif kappaonly:
        columnnames = ['Kappa Contrast']
    else:
        columnnames = ['HD Contrast', 'HR Contrast', 'HIP Contrast', 'Kappa Contrast']
elif snronly:
    if broadbandonly:
        columnnames = ['HD SNR', 'HR SNR', 'HIP SNR']
    elif kappaonly:
        columnnames = ['Kappa SNR']
    else:
        columnnames = ['HD SNR', 'HR SNR', 'HIP SNR', 'Kappa SNR']
else:
    if broadbandonly:
        columnnames = ['HD SNR', 'HD Contrast', 'HR SNR', 'HR Contrast', 'HIP SNR', 'HIP Contrast']
    elif kappaonly:
        columnnames = ['Kappa SNR', 'Kappa Contrast']
    else:
        columnnames = ['HD SNR', 'HD Contrast', 'HR SNR', 'HR Contrast', 'HIP SNR', 'HIP Contrast', 'Kappa SNR',
                       'Kappa Contrast']
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
if snronly:
    if broadbandonly:
        to_csv.to_csv('groupscores-broadband-snr.csv', index=False)
    elif kappaonly:
        to_csv.to_csv('groupscores-Kappa-snr.csv', index=False)
    else:
        to_csv.to_csv('groupscores-snr.csv', index=False)
elif contrastonly:
    if broadbandonly:
        to_csv.to_csv('groupscores-broadband-contrast.csv', index=False)
    elif kappaonly:
        to_csv.to_csv('groupscores-Kappa-contrast.csv', index=False)
    else:
        to_csv.to_csv('groupscores-contrast.csv', index=False)
else:
    if broadbandonly:
        to_csv.to_csv('groupscores-broadband.csv', index=False)
    elif kappaonly:
        to_csv.to_csv('groupscores-Kappa.csv', index=False)
    else:
        to_csv.to_csv('groupscores.csv', index=False)
