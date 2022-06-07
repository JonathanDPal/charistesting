import os
import sys
import numpy as np
from glob import glob
from valuefinder import valuefinder, paramvaluesfinder


def valuefinder2(filename, param):
    """
    Looks at a filename and discerns the KLIP parameters used to produce it. Can either find a specific KLIP
    parameter and return it in the form of a string, or it can find all KLIP parameters and return them in original
    form (int/float/bool).
    ---
    Args:
        filename (str): The name of the file we are interested in.
        param (str): Either a KLIP parameter (with the caveat that numbasis='kl' and corr_smooth='smooth'), or 'all'.
    ---
    Returns:
        If param is a specific parameter, returns a singular value for that parameter. If param is 'all',
        then returns a list of all the KLIP parameters.
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
            elif prm == 'movement' or prm == 'smooth' or prm == 'um':
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


direc = sys.argv[1]
if len(sys.argv) > 2:
    if sys.argv[2] == 'contrast':
        look_at_contrast = True
        look_at_detections = False
    elif sys.argv[2] == 'detections':
        look_at_contrast = False
        look_at_detections = True
    else:
        raise ValueError('Second keyword argument should be either "contrast" or "detections"')
else:
    look_at_contrast, look_at_detections = True, True

if look_at_contrast:
    contrast_files = glob(f'{direc}/calibrated_contrast/*.csv')
    contrast_completed = [valuefinder(file, 'all') for file in contrast_files]

if look_at_detections:
    detection_files = glob(f'{direc}/detections/*.csv')
    detection_completed = [valuefinder2(file, 'all') for file in detection_files]

olddirec = os.getcwd()
os.chdir(direc)
try:
    annuli = paramvaluesfinder('annuli')
except UnboundLocalError:
    os.chdir(olddirec)
    os.system(f'cp useful_logfile.txt {direc}/log.txt')
    os.chdir(direc)
    annuli = paramvaluesfinder('annuli')

subsections = paramvaluesfinder('subsections')
movement = paramvaluesfinder('movement')
numbasis = paramvaluesfinder('numbasis')
corr_smooth = paramvaluesfinder('corr_smooth')
highpass = paramvaluesfinder('highpass')

if direc in ['HD1160', 'HR8799', 'HIP86032', 'KB']:  # broadband
    wavelength = [float(f'{round(wvl, 2)}') for wvl in [1.15956144, 1.19969705, 1.24122187, 1.28418397, 1.32863311,
                                                        1.37462076, 1.42220017, 1.47142643, 1.52235655, 1.5750495,
                                                        1.6295663, 1.68597007, 1.74432613, 1.80470206, 1.86716776,
                                                        1.93179558, 1.99866034, 2.06783947, 2.13941309, 2.21346406,
                                                        2.29007815, 2.36934405]]
elif direc in ['Kappa']:  # k-band
    wavelength = [float(f'{round(wvl, 2)}') for wvl in [2.01513642, 2.03556326, 2.05619717, 2.07704023, 2.09809457,
                                                        2.11936233, 2.14084568, 2.1625468, 2.1844679, 2.2066112,
                                                        2.22897897, 2.25157347, 2.274397, 2.29745189, 2.32074048,
                                                        2.34426514, 2.36802826]]

if look_at_contrast:
    cc = 0
    nn = 0
    nnc = list()
    for ann in annuli:
        for sbs in subsections:
            for mov in movement:
                for spec in [None]:
                    for nb in numbasis:
                        for cs in corr_smooth:
                            for hp in highpass:
                                for wv in wavelength:
                                    if [ann, sbs, mov, spec, nb, cs, hp, wv] in contrast_completed:
                                        cc += 1
                                    else:
                                        nnc.append([ann, sbs, mov, spec, nb, cs, hp, wv])
                                        nn += 1
    print('Contrast Files:')
    print(cc)
    print(nn)

if look_at_detections:
    c = 0
    n = 0
    nc = list()
    for ann in annuli:
        for sbs in subsections:
            for mov in movement:
                for spec in [None]:
                    for nb in numbasis:
                        for cs in corr_smooth:
                            for hp in highpass:
                                if [ann, sbs, mov, spec, nb, cs, hp] in detection_completed:
                                    c += 1
                                else:
                                    nc.append([ann, sbs, mov, spec, nb, cs, hp])
                                    n += 1
    print('Detections Files:')
    print(c)
    print(n)

if look_at_detections:
    with open(f'remainingdetectionsparams.txt', 'w') as file:
        for incompleteparams in nc:
            ann, sbs, mov, spec, nb, cs, hp = incompleteparams
            file.write(f'{ann},{sbs},{mov},{spec},{nb},{cs},{hp}\n')
if look_at_contrast:
    with open(f'remainingcontrastparams.txt', 'w') as file:
        for incompleteparams2 in nnc:
            ann, sbs, mov, spec, nb, cs, hp, _ = incompleteparams2
            file.write(f'{ann},{sbs},{mov},{spec},{nb},{cs},{hp}\n')
