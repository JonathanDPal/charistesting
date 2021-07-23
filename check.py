import sys
import numpy as np
from glob import glob


def listset(alist):
    return list(set(alist))


def valuefinder(filename, param):
    """
    Looks at a filename and discerns the KLIP parameters used to produce it. Can either find a specific KLIP
    parameter and return it in the form of a string, or it can find all KLIP parameters and return them in original
    form (int/float/bool).
    """
    param = str.lower(param)
    paramlengths = {'annuli': 6, 'subsections': 11, 'movement': 8, 'spectrum': 8, 'kl': 2, 'smooth': 6, 'highpass': 8}
    if param != 'all':
        paramlength = paramlengths[param]
        startingindex = None  # will be defined soon
        for j in range(len(filename)):
            if str.lower(filename[j: j + paramlength]) == param:
                if param == 'kl':
                    startingindex = j + paramlength
                else:
                    startingindex = j - 1

        if startingindex is not None:
            if param != 'kl':
                valuelength = 0
                while startingindex >= 0 and filename[startingindex] != '_' and filename[startingindex] != '/':
                    startingindex -= 1
                    valuelength += 1
                end_index = startingindex + 1 + valuelength
                value = filename[startingindex + 1: end_index]
            else:
                end_index = startingindex + 2
                value = filename[startingindex: end_index]

        return value
    else:
        values = []
        for prm in paramlengths.keys():
            paramlength = paramlengths[prm]
            startingindex = None  # will be defined soon
            for j in range(len(filename)):
                if str.lower(filename[j: j + paramlength]) == prm:
                    if prm == 'kl':
                        startingindex = j + paramlength
                    else:
                        startingindex = j - 1

            if prm != 'kl':
                valuelength = 0
                while startingindex > 0 and filename[startingindex] != '_' and filename[startingindex] != '/':
                    startingindex -= 1
                    valuelength += 1
                end_index = startingindex + 1 + valuelength
                value = filename[startingindex + 1: end_index]
            else:
                end_index = startingindex + 2
                value = filename[startingindex: end_index]

            if prm == 'annuli' or prm == 'subsections' or prm == 'kl':
                value = int(value)
            elif prm == 'movement' or prm == 'smooth':
                value = float(value)
            elif prm == 'highpass':
                if str.lower(value) == 'true':
                    value = True
                elif str.lower(value) == 'false':
                    value = False
                else:
                    value = float(value)
            elif prm == 'spectrum':
                if str.lower(value) == 'none':
                    value = None
            values.append(value)

        return values


# direc = sys.argv[1]
direc = 'HR8799'
klip_outputs = glob(f'{direc}/klipped_cubes_Wfakes/*.fits')
completed = [valuefinder(file, 'all') for file in klip_outputs]

c = 0
n = 0
nc = list()
for ann in [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
    for sbs in [2, 4, 6]:
        for mov in [0.0, 0.5, 1.0, 1.5, 2.0]:
            for spec in [None]:
                for nb in [10, 15, 20, 25, 30, 35, 40, 50, 60]:
                    for cs in [0.0, 0.5, 1.0, 2.0]:
                        for hp in [False, 5.0, True, 15.0, 30.0]:
                            if [ann, sbs, mov, spec, nb, cs, hp] in completed:
                                c += 1
                            else:
                                nc.append([ann, sbs, mov, spec, nb, cs, hp])
                                n += 1
print(c)
print(n)