import numpy as np

annuli = [[4, 6, 8, 10, 12], [2, 3, 5, 7, 9, 11], [4, 6, 8, 10, 12], [2, 3, 5, 7, 9, 11],
          [2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12], [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
          [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
          [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
          [2, 3, 4, 5, 6, 7, 8, 9, 10, 11], [12], [12], [12]]
subsections = [[2, 4, 6], [2, 4], [2, 4], [2, 4, 6], [2, 4, 6], [2, 4, 6], [2, 4, 6], [2, 4, 6], [2, 4, 6], [2, 4, 6],
               [2, 4, 6], [6], [6], [6], [6]]
movement = [[0.0, 1.0, 2.0], [0.5, 1.5], [0.5, 1.5], [0.0, 1.0, 2.0], [0.0, 1.0, 2.0], [0.0, 1.0, 2.0],
            [0.0, 1.0, 2.0], [0.0, 1.0, 2.0], [0.5, 1.5], [0.5, 1.5], [0.5, 1.5], [0.5, 1.5], [0.5], [1.5], [1.5]]
spectrum = [[None], [None], [None], [None], [None], [None], [None], [None], [None], [None], [None], [None], [None],
            [None], [None]]
numbasis = [[10, 20, 30, 40, 50, 60], [15, 20, 25, 30, 35], [15, 20, 25, 30, 35], [10, 20, 30, 40, 50, 60],
            [15, 25, 35], [15, 25, 35], [10, 15, 20, 25, 30, 35, 40, 50, 60], [40, 10, 50, 20, 60, 30],
            [35, 40, 10, 15, 50, 20, 25, 60, 30], [35, 40, 10, 15, 50, 20, 25, 60, 30], [40, 10, 50, 60],
            [35, 15, 20, 25, 30], [35, 15, 20, 25, 30], [35, 15, 20, 25], [30]]
corr_smooth = [[0.0, 1.0, 2.0], [0.0, 0.5, 1.0], [0.0, 0.5, 1.0], [0.0, 1.0, 2.0], [0.0, 0.5, 1.0, 2.0],
               [0.0, 0.5, 1.0, 2.0], [0.0, 0.5, 1.0, 2.0], [0.5], [0.0, 1.0, 2.0, 0.5], [2.0], [0.0, 1.0, 0.5],
               [0.0, 1.0, 0.5], [0.0, 1.0, 0.5], [0.0, 1.0, 0.5], [0.0, 0.5, 1.0]]
highpass = [[False, 5.0, True, 15.0], [False, True, 30.0], [False, True, 30.0], [False, 5.0, True, 15.0],
            [False, 5.0, True, 15.0], [False, 5.0, True, 15.0], [30.0], [False, True, 5.0, 15.0], [5.0, 15.0],
            [False, True, 30.0], [False, True, 30.0], [False, True, 30.0], [False, True, 30.0], [False, True, 30.0],
            [False, True, 30.0]]

completed = []
for i in range(len(annuli)):
    ann, sbs, mov, spec, nb, cs, hp = [value[i] for value in [annuli, subsections, movement, spectrum, numbasis,
                                                                corr_smooth, highpass]]
    for a in ann:
        for s in sbs:
            for m in mov:
                for p in spec:
                    for n in nb:
                        for c in cs:
                            for h in hp:
                                completed.append([a, s, m, p, n, c, h])

fullannuli = []
fullsubsections = []
fullmovement = []
fullspectrum = []
fullnumbasis = []
fullcorrsmooth = []
fullhighpass = []
for new_full_param, fullparam in zip([fullannuli, fullsubsections, fullmovement, fullspectrum, fullnumbasis,
                                        fullcorrsmooth, fullhighpass], [annuli, subsections, movement, spectrum,
                                                                         numbasis, corr_smooth, highpass]):
    for subset in fullparam:
        new_full_param += subset


def listset(alist):
    return list(set(alist))


fullannuli, fullsubsections, fullmovement, fullspectrum, fullnumbasis, fullcorrsmooth, fullhighpass = listset(
    fullannuli), listset(fullsubsections), listset(fullmovement), listset(fullspectrum), listset(fullnumbasis), \
                                                                        listset(fullcorrsmooth), listset(fullhighpass)


fullset = []
for a in fullannuli:
    for s in fullsubsections:
        for m in fullmovement:
            for p in fullspectrum:
                for n in fullnumbasis:
                    for c in fullcorrsmooth:
                        for h in fullhighpass:
                            if not [a, s, m, p, n, c, h] in completed:
                                fullset.append([a, s, m, p, n, c, h])

with open('suggestionsHR8799.txt', 'w') as file:
    file.write(f'Annuli: {fullannuli}')
    file.write(f'\nSubsections: {fullsubsections}')
    file.write(f'\nMovement: {fullmovement}')
    file.write(f'\nSpectrum: {fullspectrum}')
    file.write(f'\nNumbasis: {fullnumbasis}')
    file.write(f'\nSmooth: {fullcorrsmooth}')
    file.write(f'\nHighpass: {fullhighpass}\n')
    for sg in fullset:
        file.write('\n'+str(sg))
