import pandas as pd
import numpy as np
import sys

snr = pd.read_csv('numericalscoring/snr_scores.csv')
contrast = pd.read_csv('numericalscoring/contrast_scores.csv')
del snr['Unnamed: 0']
del contrast['Unnamed: 0']

if len(sys.argv) > 1:
    snrweight, contrastweight = sys.argv[1], sys.argv[2]
else:
    # this should essentially allow them to be weighted equally
    snrweight, contrastweight = contrast['Score'].median(), snr['Score'].median()

overall = snr.merge(contrast, how='outer', on=['Annuli', 'Subsections', 'Movement', 'Spectrum', 'Numbasis',
                                               'Corr_Smooth', 'Highpass'], suffixes=['_snr', '_contrast'])

overall['Overall Score'] = [(s * snrweight + c * contrastweight) / (snrweight + contrastweight)
                            for s, c in zip(overall['Score_snr'], overall['Score_contrast'])]

overall.to_csv('numericalscoring/overall_scores.csv')