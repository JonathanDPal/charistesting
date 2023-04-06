import numpy as np
from glob import glob
import sys
import os
import pandas as pd

reference_contrast = [(20, 1e-5), (40, 5e-6), (60, 1e-6)]  # some values that are the standard everything is judged
# against (first value is seperation, second is standard value for that seperation)
total_injected = 18
injected_fluxes = [(20, (2e-4, 1e-4)), (40, (1e-4, 5e-5)), (60, (5e-5, 1e-5))]
detectionsfiles = glob('../detections/*.csv')


def valuefinder(filename, param):
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


reference_planets_snr = []
for reference, inj_flux in zip(reference_contrast, injected_fluxes):
    _, ref_contrast = reference
    _, fluxes = inj_flux
    for flux in fluxes:
        for _ in range(int(total_injected / (len(injected_fluxes) * len(fluxes)))):
            reference_planets_snr.append(flux * 5 / ref_contrast)
reference_score = np.sum(reference_planets_snr)

annuli, subsections, movement, numbasis, corr_smooth, highpass, score20, score40, score60, scisnr = list(), list(), \
                                                                                                    list(), list(), \
                                                                                                    list(), list(), \
                                                                                                    list(), list(), \
                                                                                                    list(), list()
for dfile in detectionsfiles:
    df = pd.read_csv(dfile)
    df1 = df[df['Injected'] == 'Science Target']
    scisnr.append(df1['SNR Value'].sum())
    df = df[df['Injected'] != 'Science Target']
    df20 = df[np.abs(df['Sep (pix)'] - 20) < 10]
    df40 = df[np.abs(df['Sep (pix)'] - 40) <= 10]
    df60 = df[np.abs(df['Sep (pix)'] - 60) < 10]

    for score, df in zip([score20, score40, score60], [df20, df40, df60]):
        cumulative_score = 0
        for snr, tf in zip(df['SNR Value'], df['Injected']):
            if type(tf) == str:
                if str.lower(tf) == 'true':
                    cumulative_score += snr
                elif str.lower(tf) == 'false':
                    cumulative_score -= snr
            else:
                if tf is True:
                    cumulative_score += snr
                elif tf is False:
                    cumulative_score -= snr
        score.append(cumulative_score / reference_score * 100)

    ann, sbs, mov, spec, nb, cs, hp = valuefinder(dfile, 'all')
    annuli.append(ann)
    subsections.append(sbs)
    movement.append(mov)
    numbasis.append(nb)
    corr_smooth.append(cs)
    highpass.append(hp)

finaldata = pd.DataFrame({'Annuli': annuli, 'Subsections': subsections, 'Movement': movement, 'Numbasis': numbasis,
                          'Corr_Smooth': corr_smooth, 'Highpass': highpass, 'Score20': score20, 'Score40': score40,
                          'Score60': score60, 'SciSNR': scisnr})

if not os.path.exists('numericalscoring'):
    os.mkdir('numericalscoring')
finaldata.to_csv('numericalscoring/snr_scores_indsep.csv', index=False)
