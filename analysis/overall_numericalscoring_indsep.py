import pandas as pd

snr = pd.read_csv('numericalscoring/snr_scores_indsep.csv')
contrast = pd.read_csv('numericalscoring/contrast_scores_indsep.csv')

overall = snr.merge(contrast, how='outer', on=['Annuli', 'Subsections', 'Movement', 'Numbasis',
                                               'Corr_Smooth', 'Highpass'])
for col in overall.columns:
    if col not in ['Annuli', 'Subsections', 'Movement', 'Numbasis', 'Corr_Smooth', 'Highpass'] and 'Score' not in col\
            and 'SciSNR' not in col:
        del overall[col]
overall.to_csv('numericalscoring/overall_scores_indsep.csv', index=False)
