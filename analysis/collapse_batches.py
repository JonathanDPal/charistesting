import pandas as pd
from glob import glob
import sys

batchfiles = glob(sys.argv[1])
df = pd.read_csv(batchfiles[0])
for bf in batchfiles:
    df = df.append(bf)
max_thold = df['Threshold'].max()
df = df[df['Threshold'] == max_thold]
df.to_csv(sys.argv[2], index=False)
