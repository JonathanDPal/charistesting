import numpy as np
import pandas as pd
from glob import glob
import sys


def valuefinder(filename, param):
    param = str.lower(param)
    paramlengths = {'annuli': 6, 'subsections': 11, 'movement': 8, 'spectrum': 8, 'smooth': 6, 'highpass': 8, 'kl': 2}
    paramlength = paramlengths[param]
    startingindex = None  # will be defined soon
    for j in range(len(filename)):
        if str.lower(filename[j: j + paramlength]) == param:
            if param == 'kl':
                startingindex = j + paramlength
            else:
                startingindex = j - 1

    if startingindex is not None:
        if param != 'kl':
            valuelength = 0
            while startingindex >= 0 and filename[startingindex] != '_' and filename[startingindex] != '/':
                startingindex -= 1
                valuelength += 1
            end_index = startingindex + 1 + valuelength
            value = filename[startingindex + 1: end_index]
        else:
            end_index = startingindex + 2
            value = filename[startingindex: end_index]

    return value


direc = sys.argv[1]
files = glob(f'{direc}/calibrated_contrast/*.csv')

ANN, SBS, MOV, NB, CS, HP = [list() for _ in range(6)]
for file in files:
    ANN.append(valuefinder(file, 'annuli'))
    SBS.append(valuefinder(file, 'subsections'))
    MOV.append(valuefinder(file, 'movement'))
    NB.append(valuefinder(file, 'kl'))
    CS.append(valuefinder(file, 'smooth'))
    HP.append(valuefinder(file, 'highpass'))

caps = [ANN, SBS, MOV, NB, CS, HP]
caps = [list(set(param)) for param in caps]

numgood, numbad = list(), list()  # how many good vs. bad vals
for file in files:
    df = pd.read_csv(file)
    good_file = True
    numgoodvals = 0
    numbadvals = 0
    for val in df['Calibrated Contrast']:
        if np.isnan(val) or val < 1e-10:
            good_file = False
            numgoodvals += 1
        else:
            numbadvals += 1
    numgood.append(numgoodvals)
    numbad.append(numbadvals)

lowers = [{val: [[], []] for val in param} for param in caps]

for i, file in enumerate(files):
    vals = [valuefinder(file, p) for p in ['annuli', 'subsections', 'movement', 'kl', 'smooth', 'highpass']]
    for val, paramdict in zip(vals, lowers):
        paramdict[val][0].append(numgood[i])
        paramdict[val][1].append(numbad[i])

ratios = [[round(float(np.sum(param_dict[key][0])) / float(np.sum(param_dict[key][1])), 4) for key in
           param_dict.keys()] for param_dict in lowers]
param_names = ['Annuli', 'Subsections', 'Movement', 'Numbasis', 'Smoothing', 'Highpass']

with open(f'{direc}_paramresults.txt', 'w') as file:
    file.write(f'{param_names[0]}:')
    for i, (c, r) in enumerate(zip(caps, ratios)):
        if i != 0:
            file.write(f'\n---\n{param_names[i]}')
        for val, ratio in zip(c, r):
            file.write(f'\n{val}: {ratio}')
