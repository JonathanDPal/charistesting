import pandas as pd
import numpy as np
import os

annuli = [4, 6, 8, 10, 12]
subsections = [2, 4, 6]
movement = [0, 1, 2]
spectrum = [None]
numbasis = [10, 20, 30, 40, 50, 60]
corr_smooth = [0, 1, 2]
highpass = [False, 5.0, True, 15.0]

filepaths = {(ann, mov): f'{ann}Annuli_{subsec}Subsections_{mov}Movement_{spec}Spectrum_{cs}Smooth_{hp}Highpass__KL'
                         '{nb}_contrast.csv' for ann in annuli for subsec in subsections for mov in movement for
             spec in spectrum for nb in numbasis for cs in corr_smooth for hp in highpass}

params = []
avgs = []
contrasts = []

for p in filepaths.keys():
    df1 = pd.read_csv(filepaths[p])
    contrast = df1['Calibrated Contrast']
    params.append(p)
    avgs.append(np.mean([c if not np.isnan(c) else 1 for c in contrast]))  # heavily punishing NaN values
    contrasts.append(contrast)

scores = [list() for _ in range(len(params))]

specific_seperations = []
for i in range(len(contrasts[0])):
    specific_seperations[i] = [c[i] for c in contrasts]

for contrast_index, ss in enumerate(specific_seperations):
    local_ss = ss
    rank = 1
    while len(local_ss) != 0:
        best_indices = [i for i, j in enumerate(local_ss) if j == np.min(local_ss)]
        for index in best_indices:
            scores[index].append(rank)
        rank += len(best_indices)
        local_ss = [local_ss[i] for i in range(len(local_ss)) if i not in best_indices]

local_avgs = avgs
rank = 1
while len(local_avgs) != 0:
    best_index_avgs = [i for i, j in enumerate(avgs) if j == np.min(local_avgs)]
    for index in best_index_avgs:
        scores[index].append(rank)
    rank += len(best_index_avgs)
    local_avgs = [avgs[i] for i in range(len(local_avgs)) if i not in best_index_avgs]
