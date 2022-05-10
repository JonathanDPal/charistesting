import numpy as np
import pandas as pd
import sys
import os

try:
    type = sys.argv[1]  # either 'snr' or 'contrast'
except IndexError:
    raise IndexError("Please specify either 'snr' or 'contrast'")

if str.lower(type) == 'snr':
    scores = pd.read_csv('analysis/numericalscoring/snr_scores.csv')
elif str.lower('contrast'):
    scores = pd.read_csv('analysis/numericalscoring/contrast_scores.csv')
else:
    raise ValueError('Check your command line argument.')


def filename(row):
    ann = row['Annuli']
    sbs = row['Subsections']
    mov = row['Movement']
    nb = row['Numbasis']
    cs = row['Corr_Smooth']
    hp = row['Highpass']
    if str.lower(hp) == 'false':
        hp = False
    else:
        hp = float(hp)
    return f'detections/{int(ann)}Annuli_{int(sbs)}Subsections_{float(mov)}Movement_NoneSpectrum_{float(cs)}Smooth' \
           f'_{hp}Highpass__KL{int(nb)}_SNR-2.csv'


targetsnr = list()
for _, r in scores.iterrows():
    df = pd.read_csv(filename(r))
    df1 = df[df['Injected'] == 'Science Target']
    if len(df1) == 0:
        targetsnr.append(0.)
    else:
        targetsnr.append(df1['SNR Value'].sum())

scores['Science SNR'] = targetsnr
direc = os.getcwd().split('/')[-2]
scores.to_csv(f'{direc}-{type}_with_target.csv', index=False)