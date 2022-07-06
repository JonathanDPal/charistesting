from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

# pandas gets mad at me if I try to do all of this in one line
df1 = pd.read_csv('numericalscoring/overall_scores.csv')
df2 = df1[df1['Score_contrast'] != -np.inf]
df3 = df2[np.isnan(df2['Score_contrast']) == False]
df = df3[np.isnan(df3['Score_snr']) == False]

sscores, cscores = df['Score_snr'], df['Score_contrast']
scisnrs = list()
analysisdirec = os.getcwd()
os.chdir('../detections')
for _, row in df.iterrows():
    ann, sbs, mov, nb, cs, hp = row['Annuli'], row['Subsections'], row['Movement'], row['Numbasis'], \
                                row['Corr_Smooth'], row['Highpass']
    filename = f'{ann}Annuli_{sbs}Subsections_{mov}Movement_NoneSpectrum_{cs}Smooth_{hp}Highpass__KL{nb}_SNR-2.csv'
    ldf = pd.read_csv(filename)
    scisnrs.append(ldf[ldf['Injected'] == 'Science Target']['SNR Value'].sum())
os.chdir(analysisdirec)

df['Sci SNR'] = scisnrs
df.to_csv('snr_contrst.csv')

sc = plt.scatter(sscores, cscores, c=scisnrs, cmap='plasma', s=4)
cbar = plt.colorbar(sc)
cbar.set_label("Sum of Science Targets' SNRs")
plt.xlabel('SNR Score')
plt.ylabel('Contrast Score')
plt.title(f'Contrast vs. SNR For {str.split(os.getcwd(), "/")[-2]}')

if len(sys.argv) == 5:
    xmin, xmax, ymin, ymax = sys.argv[1:]
    plt.xlim([int(xmin), int(xmax)])
    plt.ylim([int(ymin), int(ymax)])

plt.savefig(f'{str.split(os.getcwd(), "/")[-2]}_con_vs_snr.png')
