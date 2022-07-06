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
files = glob(f'{star}/detections/*.csv')
annuli = np.arange(12) + 1
numbasis = (np.arange(6) + 1) * 10
scores = {ann: {nb: list() for nb in numbasis} for ann in annuli}

for fle in files:
    try:
        df = pd.read_csv(fle).fillna(1)
        score = 0
        for snr, inj in zip(df['SNR Value'], df['Injected']):
            if snr < 3:
                continue
            if str(inj).lower() == 'true':
                score += snr
    except pandas.errors.EmptyDataError:
        score = 0
    ann, nb = int(valuefinder(fle, 'annuli')), int(valuefinder(fle, 'kl'))
    scores[ann][nb].append(score)

scores_median = {ann: {nb: np.median(scores[ann][nb]) for nb in numbasis} for ann in annuli}
scores_max = {ann: {nb: np.max(scores[ann][nb]) for nb in numbasis} for ann in annuli}

to_plot_median = [[scores_median[ann][nb] for nb in numbasis] for ann in annuli]
to_plot_max = [[scores_max[ann][nb] for nb in numbasis] for ann in annuli]
plotting_df_median = pd.DataFrame(to_plot_median, index=annuli, columns=numbasis)
plotting_df_max = pd.DataFrame(to_plot_max, index=annuli, columns=numbasis)

sns.heatmap(plotting_df_median, cbar=False, cmap='inferno')
plt.savefig(f'{star}/{star}_AnnNb_median_snr.png')
plt.close()

sns.heatmap(plotting_df_max, cbar=False, cmap='inferno')
plt.savefig(f'{star}/{star}_AnnNb_max_snr.png')

plotting_df_median.to_csv(f'{star}/{star}_AnnNb_median_snr.csv')
plotting_df_max.to_csv(f'{star}/{star}_AnnNb_max_snr.csv')
