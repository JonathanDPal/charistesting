import numpy as np
import sys
import pandas as pd
from math import factorial
import itertools

n = int(sys.argv[1])  # how big of parameter sets we're checking
csvfile = sys.argv[2]  # where the ranking information is coming from
numtosave = int(sys.argv[3])  # only going to save some subset of all the combinations checked
if len(sys.argv) == 5:  # can specify to only look at contrast/snr scores
    if str.lower(sys.argv[4]) == 'contrast':
        contrastonly, snronly = True, False
    elif str.lower(sys.argv[4]) == 'snr':
        contrastonly, snronly = False, True
elif len(sys.argv) == 6:  # can specify to only do Kappa/exclude Kappa
    if str.lower(sys.argv[5]) == 'contrast':
        contrastonly, snronly = True, False
    elif str.lower(sys.argv[5]) == 'snr':
        contrastonly, snronly = False, True
    else:
        contrastonly, snronly = False, False
    if str.lower(sys.argv[5]) in ['b', 'broadband']:
        broadbandonly, kappaonly = True, False
    elif str.lower(sys.argv[5]) in ['k', 'kappa']:
        broadbandonly, kappaonly = False, True
    else:
        broadbandonly, kappaonly = False, False
else:
    contrastonly, snronly, broadbandonly, kappaonly = False, False, False, False

annuli = np.arange(12) + 1
subsections = np.arange(6) + 1
movement = [float(x) for x in np.arange(7)]
corr_smooth = [float(x) for x in np.arange(4)]
highpass = ['False', '15.0', '30.0']


def choose(n, m):
    return factorial(n) / (factorial(n-m) * factorial(m))


numcombos = choose(np.prod([len(prm) for prm in [annuli, subsections, movement, corr_smooth, highpass]]), n)
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


def subdf_finder(prmcombo, pddf):
    ann, sbs, mov, cs, hp = prmcombo
    df1 = pddf[pddf['Annuli'] == ann]
    df1 = df1[df1['Subsections'] == sbs]
    df1 = df1[df1['Movement'] == mov]
    df1 = df1[df1['Corr_Smooth'] == cs]
    df1 = df1[df1['Highpass'] == hp]
    return list(df1.index)


df = pd.read_csv(csvfile)
allprmsets = [(ann, sbs, mov, cs, hp) for ann in annuli for sbs in subsections for mov in movement for cs in
              corr_smooth for hp in highpass]
allprmcombos = itertools.combinations(allprmsets, n)
data_to_save = {'Params': list(), 'Threshold': list()}
for pc in allprmcombos:
    paramgroups = [pc1 for pc1 in pc]
    data_to_save['Params'].append(str(paramgroups))
    indices = list()
    for pc1 in pc:
        indices += subdf_finder(pc1, df)
    subdf = df.iloc[indices, :]
    threshold = np.min([np.max(subdf[col]) for col in columnnames])
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
