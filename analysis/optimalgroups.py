import numpy as np
import sys
import pandas as pd
import itertools

n = int(sys.argv[1])  # how big of parameter sets we're checking
csvfile = sys.argv[2]
numrows = int(sys.argv[3])  # how many parameter sets are we using for combinations
numcombos = np.prod([numrows - k for k in range(n)]) / np.prod(np.arange(n) + 1)
print(f'{numcombos} combinations will be checked.')  # so that if it's like a trillion then I just kill the script

columnnames = ['HD SNR', 'HD Contrast', 'HR SNR', 'HR Contrast', 'HIP SNR', 'HIP Contrast', 'Kappa SNR',
               'Kappa Contrast']
df = pd.read_csv(csvfile)[:numrows]
indexcombos = itertools.combinations(np.arange(numrows), n)
data_to_save = {'Params': list(), 'Threshold': list(), 'Average': list(), 'Overall': list()}
for ic in indexcombos:
    subdf = df.iloc[list(ic), :]
    fullmin = np.min([np.min(subdf[col]) for col in columnnames])
    average = np.mean(subdf['Average'])
    paramgroups = list()
    for idx in range(len(ic)):
        paramgroups.append(tuple([df.iloc[idx, k] for k in range(7)]))
    data_to_save['Params'].append(str(paramgroups))
    data_to_save['Threshold'].append(fullmin)
    data_to_save['Average'].append(average)
    data_to_save['Overall'].append(np.mean([fullmin, fullmin, average]))

combined = pd.DataFrame(data_to_save)
to_csv = combined.sort_values('Overall', ascending=False, ignore_index=True)
to_csv.to_csv('groupscores.csv', index=False)
