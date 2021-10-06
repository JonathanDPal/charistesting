import numpy as np
from glob import glob
import sys
import pandas as pd

reference_contrast = [(20, 1e-5), (40, 5e-5), (60, 1e-6)]  # some values that are the standard everything is judged
# against (first value is seperation, second is standard value for that seperation)

contrastfiles = glob('../calibrated_contrast/*.csv')


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


reference_score = 0
for reference in reference_contrast:
    _, reference_val = reference
    reference_score += np.log10(reference_val / 5)  # stuff above is 5 sigma contrast

annuli, subsections, movement, spectrum, numbasis, corr_smooth, highpass, scores = list(), list(), list(), list(), \
                                                                                   list(), list(), list(), list()
for file in contrastfiles:
    ann, sbs, mov, spec, nb, cs, hp = valuefinder(file, 'all')
    annuli.append(ann)
    sbs.append(sbs)
    movement.append(mov)
    spectrum.append(spec)
    numbasis.append(nb)
    corr_smooth.append(cs)
    highpass.append(hp)

    df = pd.read_csv(file)
    contrast = df['Calibrated Contrast']
    seps = df['Seperation']

    score_sum = 0
    for reference in reference_contrast:
        sep, _ = reference
        closest_seperation_index = np.argmin(sep - seps)
        score_sum += (np.log10(contrast[closest_seperation_index]) / 5)  # measures 5 sigma contrast

    scores.append((score_sum / reference_score) * 100)

finaldata = pd.DataFrame({'Annuli': annuli, 'Subsections': subsections, 'Movement': movement, 'Spectrum': spectrum,
                          'Numbasis': numbasis, 'Corr_Smooth': corr_smooth, 'Highpass': highpass, 'Score': scores})
sorted_by_score = finaldata.sort_values(by='Score', ascending=False)
if not os.path.exists('numericalscoring'):
    os.mkdir('numericalscoring')
sorted_by_score.to_csv('numericalscoring/contrast_scores.csv')
