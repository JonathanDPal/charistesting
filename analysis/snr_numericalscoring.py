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
    ref_contrast /= 5  # stuff above is 5 sigma contrast
    _, fluxes = inj_flux
    for flux in fluxes:
        for _ in range(int(total_injected / (len(injected_fluxes) * len(fluxes)))):
            reference_planets_snr.append(flux / ref_contrast)
snr_values = np.arange(start=2, stop=np.ceil(np.max(reference_planets_snr)), step=1)
ref_individual_scores = list()
for snr in snr_values:
    tp = len([ref_snr for ref_snr in reference_planets_snr if ref_snr >= snr])
    ind_score = snr * (tp / total_injected)
    ref_individual_scores.append(ind_score)
reference_score = np.mean(ref_individual_scores)

annuli, subsections, movement, spectrum, numbasis, corr_smooth, highpass, score = list(), list(), list(), list(), \
                                                                                  list(), list(), list(), list()
for dfile in detectionsfiles:
    df = pd.read_csv(dfile)
    df = df[df['Injected'] != 'Science Target']  # ignoring science targets for scoring
    try:
        snrvals = np.arange(start=2, stop=np.ceil(df['SNR Value'].max()), step=1)
    except ValueError:  # this means that the detections file is empty
        continue

    individual_scores = list()
    for snr in snrvals:
        df1 = df[df['SNR Value'] >= snr]
        inj = list(df1['Injected'])
        tp = 0
        fp = 0
        for tf in inj:
            if type(tf) == str:
                if str.lower(tf) == 'true':
                    tp += 1
                elif str.lower(tf) == 'false':
                    fp += 1
            else:
                if tf is True:
                    tp += 1
                elif tf is False:
                    fp += 1
        if tp == 0 and fp == 0:
            continue
        elif tp == 0 and fp > 0:
            ind_score = -1 * snr * fp / total_injected
            individual_scores.append(ind_score)
        else:
            ind_score = snr * (tp / total_injected) * ((tp - fp) / (tp + fp))
            individual_scores.append(ind_score)

    cumulative_score = np.mean(individual_scores)
    score.append(cumulative_score / reference_score * 100)

    ann, sbs, mov, spec, nb, cs, hp = valuefinder(dfile, 'all')
    annuli.append(ann)
    subsections.append(sbs)
    movement.append(mov)
    spectrum.append(spec)
    numbasis.append(nb)
    corr_smooth.append(cs)
    highpass.append(hp)

finaldata = pd.DataFrame({'Annuli': annuli, 'Subsections': subsections, 'Movement': movement, 'Spectrum': spectrum,
                          'Numbasis': numbasis, 'Corr_Smooth': corr_smooth, 'Highpass': highpass, 'Score': score})
sorted_by_score = finaldata.sort_values(by='Score', ascending=False)
del sorted_by_score['Unnamed: 0']

if not os.path.exists('numericalscoring'):
    os.mkdir('numericalscoring')
sorted_by_score.to_csv('numericalscoring/snr_scores.csv')
