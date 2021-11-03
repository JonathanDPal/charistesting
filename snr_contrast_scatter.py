import astropy.time.utils
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os

# pandas gets mad at me if I try to do all of this in one line
df1 = pd.read_csv('numericalscoring/overall_scores.csv')
df2 = df1[df1['Score_contrast'] != -np.inf]
df3 = df2[np.isnan(df2['Score_contrast']) == False]
df = df3[np.isnan(df3['Score_snr']) == False]

sscores, cscores = df['Score_snr'], df['Score_contrast']
plt.scatter(sscores, cscores)
plt.xlabel('SNR Score')
plt.ylabel('Contrast Score')
plt.title(f'Contrast vs. SNR For {str.split(os.getcwd(), "/")[-2]}')
plt.savefig(f'{str.split(os.getcwd(), "/")[-2]}_con_vs_snr.png')
