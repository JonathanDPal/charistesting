import pandas as pd
import numpy as np
import sys

df = pd.read_csv(f'{sys.argv[1]}/analysis/numericalscoring/overall_scores.csv')
df[f'Score_{sys.argv[1]}'] = df['Score']
del df['Unnamed: 0']
del df['Score']

for dataset in sys.argv[2:]:
    newdf = pd.read_csv(f'{dataset}/analysis/numericalscoring/overall_scores.csv')
    newdf[f'Score_{dataset}'] = newdf['Score']
    del newdf['Unnamed: 0']
    del newdf['Score_snr']
    del newdf['Score_contrast']
    del newdf['Score']
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
df.to_csv('compiled_overall_numericalscores.csv')
