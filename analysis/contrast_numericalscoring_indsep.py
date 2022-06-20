import numpy as np
from glob import glob
import sys
import pandas as pd
import pandas.errors
import os
import warnings

warnings.filterwarnings("error", category=RuntimeWarning)  # error if dividing by zero or taking log of negative number

reference_contrast = [(13, 1.2e-5), (30, 2e-6), (50, 1e-6)]  # some values that are the standard everything is judged
# against (first value is seperation, second is standard value for that seperation). standard values taken from an
# SPIE proceeding paper

if len(sys.argv) == 2:
    if sys.argv[1] == 'all':
        wavelength = ''
    else:
        wavelength = sys.argv[1]  # in microns
else:
    wavelength = 1.63  # what we've been looking at on broadband
contrastfiles = glob(f'../calibrated_contrast/*{wavelength}um*.csv')
try:
    assert len(contrastfiles) != 0
except AssertionError:
    raise ValueError('The default wavelength (1.63 um) is not a wavelength of the dataset that is being analyzed '
                     'here. Please specify the wavelength on the command line, i.e. say "python contrast_numerical '
                     'scoring.py {insert wavelength here}"')


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


annuli, subsections, movement, numbasis, corr_smooth, highpass, scores13, scores30, scores50 = list(), list(), \
                                                                                               list(), list(), \
                                                                                               list(), list(), \
                                                                                               list(), list(), list()
for cfile in contrastfiles:
    try:
        df = pd.read_csv(cfile)
    except pandas.errors.EmptyDataError:
        print(f'{cfile} had no columns to parse. Will not be included in final CSV')
        continue
    contrast = df['Calibrated Contrast']
    seps = df['Seperation']

    for score_lst, reference in zip([scores13, scores30, scores50], reference_contrast):
        score_sum = 0
        sep, reference_val = reference
        closest_seperation_index = np.argmin(sep - seps)
        if contrast[closest_seperation_index] == -np.inf:
            score_sum += -np.inf
        else:
            try:
                ratio = (np.log10(contrast[closest_seperation_index])) / (np.log10(reference_val))
                score_sum += ratio
            except RuntimeWarning:  # this has happened when a negctive number is the value for contrast
                score_sum += -np.inf
        score_lst.append(score_sum / len(reference_contrast) * 100)

    ann, sbs, mov, spec, nb, cs, hp = valuefinder(cfile, 'all')
    annuli.append(ann)
    subsections.append(sbs)
    movement.append(mov)
    numbasis.append(nb)
    corr_smooth.append(cs)
    highpass.append(hp)

finaldata = pd.DataFrame({'Annuli': annuli, 'Subsections': subsections, 'Movement': movement, 'Numbasis': numbasis,
                          'Corr_Smooth': corr_smooth, 'Highpass': highpass, 'Score13': scores13, 'Score30': scores30,
                          'Score50': scores50})

if len(sys.argv) > 1 and sys.argv[1] == 'all':  # need to collapse all wavelength scores into a single contrast score
    d13, d30, d50 = dict(), dict(), dict()
    for _, row in finaldata.iterrows():
        idx = tuple(row[:-3])
        if idx in d13.keys():
            d13[idx].append(row[-1])
        else:
            d13[idx] = [row[-1]]
        if idx in d30.keys():
            d30[idx].append(row[-2])
        else:
            d30[idx] = [row[-2]]
        if idx in d50.keys():
            d50[idx].append(row[-3])
        else:
            d50[idx] = [row[-3]]
    d13 = {key: np.mean(d13[key]) for key in d13.keys()}
    d30 = {key: np.mean(d30[key]) for key in d30.keys()}
    d50 = {key: np.mean(d50[key]) for key in d50.keys()}
    annuli, subsections, movement, numbasis, corr_smooth, highpass, sco13, sco30, sco50 = list(), list(), list(), \
                                                                                          list(), list(), list(), \
                                                                                          list(), list(), list()
    for key in d13.keys():
        annuli.append(key[0])
        subsections.append(key[1])
        movement.append(key[2])
        numbasis.append(key[3])
        corr_smooth.append(key[4])
        highpass.append(key[5])
        sco13.append(d13[key])
        sco30.append(d30[key])
        sco50.append(d50[key])
    finaldata = pd.DataFrame({'Annuli': annuli, 'Subsections': subsections, 'Movement': movement, 'Numbasis':
                               numbasis, 'Corr_Smooth': corr_smooth, 'Highpass': highpass, 'Score13': sco13,
                              'Score30': sco30, 'Score50': d50})
sorted_by_score = finaldata.sort_values(by='Score', ascending=False, ignore_index=True)

if not os.path.exists('numericalscoring'):
    os.mkdir('numericalscoring')
sorted_by_score.to_csv('numericalscoring/contrast_scores_indsep.csv', index=False)
