import sys
import pandas as pd
from glob import glob

direc = sys.argv[1]
files = glob(f'{direc}/*.xlsx')

good = 0
bad = 0

for file in files:
	badfile = False
	df = pd.read_excel(file)
	for col in df.columns:
		if 'Retrieved Flux' in col:
			for elm in df[col]:
				if elm > 0:
					good += 1
				else:
					bad += 1
					badfile = True

print(f'positive retrieved fluxes in {round(float(good) / (float(good) + float(bad)) * 100, 2)}% of cases')
