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

contrastfiles = glob(f'../calibrated_contrast/*.csv')
assert len(contrastfiles) in [798336, 616896]
object_filtnames = {'HD1160': 'Broadband', 'HR8799': 'Broadband', 'KB': 'Broadband', 'Kappa': 'K'}
object_name = os.getcwd().split('/')[-2]


def find_bin_weights(filt):
    """
    Taylor's code.
    """
    filt = filt.lower()
    # Spectral resolution, R, for CHARIS lowres (broadband) or hires (J, H, K) modes
    CHARIS_spec_res = {'lowres': 30, 'hires': 100}

    # Sets ends of wavelength bands, in um, for each band; uses same values as buildcal code, which are used to set the
    #   values in the extractcube cubes that use those cals
    CHARIS_filter_ends = {'j': [1.155, 1.340], 'h': [1.470, 1.800], 'k': [2.005, 2.380], 'broadband': [1.140, 2.410]}

    # Calculating edges and midpoints of wavelength bins #

    # First get R from dictionary based on filter
    if filt in ['j', 'h', 'k']:
        R = CHARIS_spec_res['hires']
    elif filt == 'broadband':
        R = CHARIS_spec_res['lowres']
    else:
        raise ValueError("Filter {0} not recognized. Please enter 'J', 'H', 'K', or 'Broadband'".format(filt))

    # Calculate wavelength bin end and midpoints like buildcal does
    Nspec = int(np.log(CHARIS_filter_ends[filt][1] / CHARIS_filter_ends[filt][0]) * R + 1.5)
    loglam_endpts = np.linspace(np.log(CHARIS_filter_ends[filt][0]), np.log(CHARIS_filter_ends[filt][1]), Nspec)
    lam_endpts = np.exp(loglam_endpts)

    # Calculating the wavelength-averaged contrast #

    # First, the normalized weight for each contrast is the fraction of the total wavelength range covered
    #   that is contained in that wavelength bin (note: the 'divide by the total' part of the average is
    #   already contained in the weights. If each wavelength bin was the same width, the weight would be
    #   1 / Number_of_bins )

    bin_widths = (lam_endpts[1:] - lam_endpts[:-1])
    bin_weights = bin_widths / (CHARIS_filter_ends[filt][1] - CHARIS_filter_ends[filt][0])
    return bin_weights
    # # Then just multiply each contrast by the corresponding weight and sum together
    # mean_contrast = np.sum(bin_weights * contrast_spec)


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


annuli, subsections, movement, numbasis, corr_smooth, highpass, wavelength, \
con13, con30, con50 = list(), list(), list(), list(), list(), list(), list(), list(), list(), list()

for cfile in contrastfiles:
    try:
        df = pd.read_csv(cfile)
    except pandas.errors.EmptyDataError:
        print(f'{cfile} had no columns to parse. Will not be included in final CSV')
        continue
    contrast = df['Calibrated Contrast']
    seps = df['Seperation']

    for con_lst, reference in zip([con13, con30, con50], reference_contrast):
        sep, reference_val = reference
        closest_seperation_index = np.argmin(np.abs(sep - seps))
        con_lst.append(contrast[closest_seperation_index])

    ann, sbs, mov, spec, nb, cs, hp = valuefinder(cfile, 'all')
    annuli.append(ann)
    subsections.append(sbs)
    movement.append(mov)
    numbasis.append(nb)
    corr_smooth.append(cs)
    highpass.append(hp)
    wavelength.append(float(cfile.split('_')[-2][:-2]))


finaldata = pd.DataFrame({'Annuli': annuli, 'Subsections': subsections, 'Movement': movement, 'Numbasis': numbasis,
                          'Corr_Smooth': corr_smooth, 'Highpass': highpass, 'Wavelength': wavelength,
                          'Con13': con13, 'Con30': con30, 'Con50': con50}).sort_values(['Annuli', 'Subsections',
                                                                                        'Movement', 'Numbasis',
                                                                                        'Corr_Smooth', 'Highpass',
                                                                                        'Wavelength'])

if len(sys.argv) > 1 and sys.argv[1] == 'all':  # need to collapse all wavelength scores into a single contrast score
    bin_weights = find_bin_weights(object_filtnames[object_name])
    d13, d30, d50 = dict(), dict(), dict()
    for _, row in finaldata.iterrows():
        idx = tuple(row[:-3])
        if idx in d13.keys():
            d13[idx].append(row['Con13'])
        else:
            d13[idx] = [row['Con13']]
        if idx in d30.keys():
            d30[idx].append(row['Con30'])
        else:
            d30[idx] = [row['Con30']]
        if idx in d50.keys():
            d50[idx].append(row['Con50'])
        else:
            d50[idx] = [row['Con50']]
    d13 = {key: -np.log(np.sum(np.array(d13[key]) * bin_weights)) if -np.inf not in np.array(d13[key]) else -np.inf
           for key in d13.keys()}
    d30 = {key: -np.log(np.sum(np.array(d30[key]) * bin_weights)) if -np.inf not in np.array(d30[key]) else -np.inf
           for key in d13.keys()}
    d50 = {key: -np.log(np.sum(np.array(d50[key]) * bin_weights)) if -np.inf not in np.array(d50[key]) else -np.inf
           for key in d13.keys()}
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
                              'Score30': sco30, 'Score50': sco50})

if not os.path.exists('numericalscoring'):
    os.mkdir('numericalscoring')
finaldata.to_csv('numericalscoring/contrast_scores_indsep.csv', index=False)
