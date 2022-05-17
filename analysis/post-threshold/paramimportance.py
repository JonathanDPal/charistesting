import pandas as pd
import numpy as np
import itertools


def findem(df, pms):
    p1, p2, p3, p4, p5 = pms
    k1, v1 = p1
    df1 = df[df[k1] == v1]
    for prm in [p2, p3, p4, p5]:
        key, value = prm
        df1 = df1[df1[key] == value]
    return df1


df = pd.read_csv('~/Desktop/soloscoresall.csv')
annuli = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
subsections = [1, 2, 3, 4, 5, 6]
movement = [0, 1, 2, 3, 4, 5, 6]
numbasis = [10, 20, 30, 40, 50, 60]
corr_smooth = [0, 1, 2, 3]
highpass = ['False', '15.0', '30.0']

anns, sbss, movs, nbs, css, tshs = list(), list(), list(), list(), list(), list()
for ann in annuli:
    for sbs in subsections:
        for mov in movement:
            for nb in numbasis:
                for cs in corr_smooth:
                    prms = (('Annuli', ann), ('Subsections', sbs), ('Movement', mov), ('Numbasis', nb),
                            ('Corr_Smooth', cs))
                    df2 = findem(df=df, pms=prms)
                    tshs.append(df2['Threshold'].max() - df2['Threshold'].std())
                    anns.append(ann)
                    sbss.append(sbs)
                    movs.append(mov)
                    nbs.append(nb)
                    css.append(cs)

to_save = pd.DataFrame({'Annuli': anns, 'Subsections': sbss, 'Movement': movs, 'Numbasis': nbs, 'Corr_Smooth': css,
                        'stdev': tshs})
to_save.to_csv('tshdiffs.csv', index=False)
