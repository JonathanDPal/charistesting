import numpy as np
import pandas as pd
import sys

starsystems = sys.argv[1:]
nums = [13, 20, 30, 40, 50, 60]
df = pd.read_csv(f'{starsystems[0]}/analysis/numericalscoring/overall_scores_indsep.csv')
for num in nums:
    df[f'{starsystems[0]}_{num}'] = df[f'Score{num}'].rank(pct=True)
    del df[f'Score{num}']
df[f'{starsystems[0]}_SciSNR'] = df['SciSNR'].rank(pct=True)
del df['SciSNR']

for star in starsystems[1:]:
    newdf = pd.read_csv(f'{star}/analysis/numericalscoring/overall_scores_indsep.csv')
    for num in nums:
        newdf[f'{star}_{num}'] = newdf[f'Score{num}'].rank(pct=True)
        del newdf[f'Score{num}']
    newdf[f'{star}_SciSNR'] = newdf['SciSNR'].rank(pct=True)
    del newdf['SciSNR']
    for col in newdf.columns:
        if col not in ['Annuli', 'Subsections', 'Movement', 'Numbasis', 'Corr_Smooth', 'Highpass'] + \
                        [f'{star}_{num}' for num in nums] + [f'{star}_SciSNR']:
            del newdf[col]
    df = df.merge(newdf, how='outer', on=['Annuli', 'Subsections', 'Movement', 'Numbasis', 'Corr_Smooth', 'Highpass'])

avgs = list()
for _, row in df.iterrows():
    scores = [row[col] for col in df.columns if np.sum([f'{num}' in col for num in nums]) > 0]
    scores += [row[col] for col in [f'{star}_SciSNR' for star in starsystems]]
    avg = np.mean(scores)
    avgs.append(avg)
df['Average'] = avgs

df = df.sort_values(by='Average', ascending=False, ignore_index=True)
df.to_csv('paramrankings/paramrankings_indsep.csv', index=False)
