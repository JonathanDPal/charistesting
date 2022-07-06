import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import sys
from glob import glob
from valuefinder import valuefinder
import pandas.errors

star = sys.argv[1]
con_files = glob(f'{star}/calibrated_contrast/*.csv')

df0 = pd.read_csv(con_files[0])
seps = np.array(df0['Seperation'][2:])

movlsts = [list() for _ in range(7)]

for fle in con_files:
    try:
        df1 = pd.read_csv(fle).fillna(1)
        contrast = np.array(df1['Calibrated Contrast'][2:])
    except pandas.errors.EmptyDataError:
        contrast = np.array([1 for _ in range(len(seps))])
    mov = int(float(valuefinder(fle, 'movement')))
    assert mov in [0, 1, 2, 3, 4, 5, 6]
    if len(contrast) == len(seps):
        movlsts[mov].append(contrast)

movlsts = [[np.median([ml[i] for ml in movlst]) for i in range(len(seps))] for movlst in movlsts]
df = pd.DataFrame(data=movlsts).transpose()
df.columns = [f'Mov {num}' for num in range(7)]
df.index = seps
df.to_csv(f'{star}/{star}_mov_contrast.csv')

for mov, movlst in enumerate(movlsts):
    plt.plot(seps, movlst, label=mov)

plt.xlabel('Seperation')
plt.ylabel('Median Calibrated Contrast')
plt.yscale('log')
plt.legend()
plt.savefig(f'{star}/mov_contrast.png')
