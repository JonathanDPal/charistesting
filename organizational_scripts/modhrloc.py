import pandas as pd
import sys
from math import sqrt
from glob import glob

files = glob(sys.argv[1])
for f in files:
    df = pd.read_csv(f)
    for _, row in df.iterrows():
        if str(row['Injected']).lower() == 'science target' and \
                sqrt((row['row'] - 137) ** 2 + (row['col'] - 157) ** 2) < 4:
            row['Injected'] = 'False'
        if sqrt((row['row'] - 107) ** 2 + (row['col'] - 123) ** 2) < 4:
            row['Science Target'] = 'True'
    df.to_csv(f, index=False)
