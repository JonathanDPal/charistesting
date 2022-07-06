import matplotlib.pyplot as plt
import pandas as pd
import pandas.errors
import numpy as np
from matplotlib import pyplot
import sys
from glob import glob
from valuefinder import valuefinder
import seaborn as sns

star = sys.argv[1]
starname = sys.argv[2]
df = pd.read_csv(f'{star}/analysis/numericalscoring/overall_scores.csv')
annuli = np.arange(12) + 1
numbasis = (np.arange(6) + 1) * 10
scores = {ann: {nb: list() for nb in numbasis} for ann in annuli}

for _, row in df.iterrows():
    ann, nb = row['Annuli'], row['Numbasis']
    scores[ann][nb].append(row['Score_contrast'])

scores_median = {ann: {nb: np.nanmedian(scores[ann][nb]) for nb in numbasis} for ann in annuli}
scores_max = {ann: {nb: np.nanmax(scores[ann][nb]) for nb in numbasis} for ann in annuli}

to_plot_median = [[scores_median[ann][nb] for nb in numbasis] for ann in annuli]
to_plot_max = [[scores_max[ann][nb] for nb in numbasis] for ann in annuli]
plotting_df_median = pd.DataFrame(to_plot_median, index=annuli, columns=numbasis)
plotting_df_max = pd.DataFrame(to_plot_max, index=annuli, columns=numbasis)

sns.heatmap(plotting_df_median, cbar=True, cmap='inferno', cbar_kws={'label': 'percentage of reference score'})
plt.xlabel('Number of KL Modes')
plt.ylabel('Number of Annuli')
plt.title(f'Median Contrast Score on {starname}')
plt.savefig(f'{star}/{star}_AnnNb_median.png')
plt.close()

sns.heatmap(plotting_df_max, cbar=True, cmap='inferno', cbar_kws={'label': 'percentage of reference score'})
plt.xlabel('Number of KL Modes')
plt.ylabel('Number of Annuli')
plt.title(f'Top Contrast Score on {starname}')

plt.savefig(f'{star}/{star}_AnnNb_max.png')

plotting_df_median.to_csv(f'{star}/{star}_AnnNb_median.csv')
plotting_df_max.to_csv(f'{star}/{star}_AnnNb_max.csv')
