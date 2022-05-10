import numpy as np
import pandas as pd

snrscores = pd.read_csv('analysis/numericalscoring/snr_scores.csv')


def filename(row):
    ann = row['Annuli']
    sbs = row['Subsections']
    mov = row['Movement']
    nb = row['Numbasis']
    cs = row['Corr_Smooth']
    hp = row['Highpass']
    return f'detections/{ann}Annuli_{sbs}Subsections_{mov}Movement_NoneSpectrum_{cs}Smooth_{hp}Highpass__KL' \
           f'{nb}_SNR-2.csv'


targetsnr = list()
for r in snrscores:
    df = pd.read_csv(filename(r))
    df1 = df[df['Injected'] == 'Science Target']
    if len(df1) == 0:
        targetsnr.append(0.)
    else:
        targetsnr.append(df1['SNR Value'].sum())

snrscores['Science SNR'] = targetsnr
snrscores.to_csv('snr_with_target.csv')