import pandas as pd
import numpy as np

annuli = [6, 7, 8, 9, 10, 11, 12]
subsections = [4]
movement = [0,1,2,3,4,5]
spectrum = [None]
numbasis = [20]
corr_smooth = [1]
highpass = [True]

filepaths = {(ann, mov) : '{0}Annuli_{1}Subsections_{2}Movement_{3}Spectrum_{4}Smooth_{5}Highpass__KL'
                          '{6}_contrast.csv'.format(ann, subsec, mov, spec, cs, hp, nb) for ann in annuli for subsec
             in subsections for mov in movement for spec in spectrum for nb in numbasis for cs in corr_smooth for hp
             in highpass}

params = []
sums = []
avgs = []
contrasts = []

for p in filepaths.keys():
    df1 = pd.read_csv(filepaths[p])
    contrast = df1['Calibrated Contrast']
    params.append(p)
    sums.append(np.sum(contrast))
    avgs.append(np.mean(contrast))
    contrasts.append(contrast)

best_contrast_sum = np.min(sums)
best_index_sum = [i for i, j in enumerate(sums) if j == best_contrast_sum]
for i in best_index_sum:
    print("By minimum of sum:")
    print(params[i])
    print(contrasts[i])
    print("")

best_contrast_avg = np.min(avgs)
best_index_avg = [i for i, j in enumerate(avgs) if j == best_contrast_avg]
for i in best_index_avg:
    print("By minimum of avg:")
    print(params[i])
    print(contrasts[i])


