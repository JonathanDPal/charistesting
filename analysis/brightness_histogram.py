from matplotlib import pyplot as plt
import pandas as pd
import sys
from glob import glob
import numpy as np

bright, dim = list(), list()
for direc in sys.argv[1:]:
    detectionsfiles = glob(f'{direc}/detections/*.csv')
    for f in detectionsfiles:
        df = pd.read_csv(f)
        df1 = df[df['Injected'] == 'True']
        bgt = list(df1[np.min(np.array((19, 139, 259)) - (df1['PA'])) < 3]['SNR Value'])
        dm = list(df1[np.min(np.array((79, 199, 319)) - df1['PA'] - 40) < 3]['SNR Value'])

        for lst in (bgt, dm):
            if len(lst) < 9:
                for _ in range(9 - len(lst)):
                    lst.append(0)

        bright.append(bgt)
        dim.append(dm)

bright = np.array(bright, dtype=object).flatten()
dim = np.array(dim, dtype=object).flatten()

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
ax1.hist(bright)
ax1.set_title('Brighter Tier')
ax2.hist(dim)
ax2.set_title('Dimmer Tier')
plt.savefig('tierVsnr_histogram.png')
