from matplotlib import pyplot as plt
import pandas as pd
import sys
from glob import glob
import numpy as np

sep20, sep40, sep60 = list(), list(), list()
for direc in sys.argv[1:]:
    detectionsfiles = glob(f'{direc}/detections/*.csv')
    for f in detectionsfiles:
        df = pd.read_csv(f)
        df1 = df[df['Injected'] == 'True']
        at20 = list(df1[(df1['Sep (pix)'] - 20) < 1]['SNR Value'])
        at40 = list(df1[(df1['Sep (pix)'] - 40) < 1]['SNR Value'])
        at60 = list(df1[(df1['Sep (pix)'] - 40) < 1]['SNR Value'])

        for lst in (at20, at40, at60):
            if len(lst) < 6:
                for _ in range(6 - len(lst)):
                    lst.append(0)

        sep20.append(at20)
        sep40.append(at40)
        sep60.append(at60)

sep20 = np.array(sep20, dtype=object).flatten()
sep40 = np.array(sep40, dtype=object).flatten()
sep60 = np.array(sep60, dtype=object).flatten()

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
ax1.hist(sep20)
ax1.set_title('Sep 20')
ax2.hist(sep40)
ax2.set_title('Sep 40')
ax3.hist(sep60)
ax3.set_title('Sep 60')
plt.savefig('sepVsnr_histogram.png')
