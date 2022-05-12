import numpy as np
import pandas as pd

starsystems = sys.argv[1:]
df = pd.read_csv(f'{starsystems[0]}/analysis/numericalscoring/overall_scores.csv')
df[f'{starsystems[0]}_snr'] = df['Score_snr'].rank(pct=True)
df[f'{starsystems[0]}_contrast'] = df['Score_contrast'].rank(pct=True)
del df['Score_snr']
del df['Score_contrast']
del df['Overall Score']

for star in starsystems[1:]:
    newdf = pd.read_csv(f'{star}/analysis/numericalscoring/overall_scores.csv')
    newdf[f'{star}_snr'] = newdf['Score_snr'].rank(pct=True)
    newdf[f'{star}_contrast'] = newdf['Score_contrast'].rank(pct=True)
    del newdf['Score_snr']
    del newdf['Score_contrast']
    del newdf['Overall Score']
    df = df.merge(newdf, how='outer', on=['Annuli', 'Subsections', 'Movement', 'Numbasis', 'Corr_Smooth', 'Highpass'])

avgs = list()
for _, row in df.iterrows():
    scores = [row[f'{star}_snr'] for star in starsystems] + [row[f'{star}_contrast'] for star in starsystems]
    avg = np.mean(scores)
    avgs.append(avg)
df['Average'] = avgs

df = df.sort_values(by='Average', ascending=False, ignore_index=True)
df.to_csv('paramrankings.csv', index=False)
