import pandas as pd
import numpy as np
import itertools
import sys

df = pd.read_csv('paramrankings_indsep.csv')
annuli = np.arange(12) + 1
subsections = np.arange(6) + 1
movement = [float(x) for x in np.arange(6)]
corr_smooth = [float(x) for x in np.arange(4)]
highpass = ['False', '15.0', '30.0']
newdf = pd.DataFrame(columns=[col for col in df.columns if col != 'Numbasis'])
numnb_touse = int(sys.argv[1])

stars, nums = ['HD1160', 'HR8799', 'KB', 'Kappa'], [13, 20, 30, 40, 50, 60, 'SciSNR']
columnnames = [f'{star}_{num}' for star in stars for num in nums]
pcols = ['Annuli', 'Subsections', 'Movement', 'Corr_Smooth', 'Highpass']

for ann in annuli:
    for sbs in subsections:
        for mov in movement:
            for cs in corr_smooth:
                for hp in highpass:
                    subdf = df[df['Annuli'] == ann]
                    subdf = subdf[subdf['Subsections'] == sbs]
                    subdf = subdf[subdf['Movement'] == mov]
                    subdf = subdf[subdf['Corr_Smooth'] == cs]
                    subdf = subdf[subdf['Highpass'] == hp]
                    assert len(subdf) == 6
                    cbs = [cbo for cbo in itertools.combinations(np.arange(6), numnb_touse)]
                    tholds = [np.min([subdf.iloc[list(cbo), :][col].max() for col in columnnames]) for cbo in cbs]
                    ssdf = subdf.iloc[list(cbs[np.argmax(tholds)]), :]
                    assert len(ssdf) == numnb_touse
                    newrow = [ssdf[col][ssdf.index[0]] for col in pcols]
                    maxes = [ssdf[col].max() for col in columnnames]
                    newrow += maxes
                    newrow += [np.mean(maxes)]
                    newdf.loc[len(newdf.index)] = newrow

newdf = newdf.sort_values('Average', ascending=False)
newdf.to_csv(f'paramrankings/paramrankings_indsep_nb{numnb_touse}.csv', index=False)
