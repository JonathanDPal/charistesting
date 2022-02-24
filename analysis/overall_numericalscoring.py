import pandas as pd
import numpy as np
import sys

snr = pd.read_csv('numericalscoring/snr_scores.csv')
contrast = pd.read_csv('numericalscoring/contrast_scores.csv')

if len(sys.argv) > 1:
    snrweight, contrastweight = float(sys.argv[1]), float(sys.argv[2])
else:
    # this should essentially allow them to be weighted (roughly) equally
    snrweight, contrastweight = contrast['Score'].median(), snr['Score'].median()

overall = snr.merge(contrast, how='outer', on=['Annuli', 'Subsections', 'Movement', 'Spectrum', 'Numbasis',
                                               'Corr_Smooth', 'Highpass'], suffixes=['_snr', '_contrast'])

overall['Overall Score'] = [(s * snrweight + c * contrastweight) / (snrweight + contrastweight)
                            for s, c in zip(overall['Score_snr'], overall['Score_contrast'])]

sorted_by_score = overall.sort_values(by='Overall Score', ascending=False, ignore_index=True)
sorted_by_score.to_csv('numericalscoring/overall_scores.csv', index=False)