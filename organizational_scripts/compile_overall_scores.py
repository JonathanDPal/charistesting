import pandas as pd
import numpy as np
import sys

df = pd.read_csv(f'{sys.argv[1]}/analysis/numericalscoring/overall_scores.csv')
df[f'Score_{sys.argv[1]}'] = df['Overall Score']
del df['Score_snr']
del df['Score_contrast']
del df['Overall Score']

for dataset in sys.argv[2:]:
    newdf = pd.read_csv(f'{dataset}/analysis/numericalscoring/overall_scores.csv')
    newdf[f'Score_{dataset}'] = newdf['Overall Score']
    del newdf['Score_snr']
    del newdf['Score_contrast']
    del newdf['Overall Score']
    df = df.merge(newdf, how='outer', on=['Annuli', 'Subsections', 'Movement', 'Spectrum', 'Numbasis',
                                               'Corr_Smooth', 'Highpass'])

totalscore = []
for _, row in df.iterrows():
    total = 0
    for col in df.columns:
        if 'Score' in col:
            total += row[col]
    totalscore.append(total)

df['Compiled Score'] = totalscore

sorted_by_score = df.sort_values(by='Compiled Score', ascending=False)
del sorted_by_score['Unnamed: 0']
sorted_by_score.to_csv('compiled_overall_numericalscores.csv')
