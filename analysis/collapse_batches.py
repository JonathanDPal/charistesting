import pandas as pd
from glob import glob
import sys
import os

batchfiles = glob(sys.argv[1])
df = pd.read_csv(batchfiles[0])
for bf in batchfiles[1:]:
    df = df.append(bf)
max_thold = df['Threshold'].max()
df = df[df['Threshold'] == max_thold]
df.to_csv(sys.argv[2], index=False)
for bf in batchfiles:
    os.remove(bf)
