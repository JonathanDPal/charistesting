import pandas as pd
import sys

snr_old = pd.read_csv('numericalscoring/snr_scores.csv')
snr_new = pd.read_csv('numericalscoring/snr_scores_new.csv')

merged = snr_old.merge(snr_new, how='outer', on=['Annuli', 'Subsections', 'Movement', 'Spectrum', 'Numbasis',
                                                 'Corr_Smooth', 'Highpass'], suffixes=['_old', '_new'])

merged.to_csv('numericalscoring/merged_snr.csv')
merged.to_csv(f'../../to_mac/merged_snr_{sys.argv[1]}.csv')