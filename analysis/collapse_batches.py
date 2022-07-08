import pandas as pd
from glob import glob
import sys
import os

batchfiles = glob(sys.argv[1])
df = pd.read_csv(batchfiles[0])
for bf in batchfiles[1:]:
    df = pd.concat([df, pd.read_csv(bf)], ignore_index=True)
max_thold = df['Threshold'].max()
df = df[df['Threshold'] == max_thold]
df.to_csv(sys.argv[2], index=False)