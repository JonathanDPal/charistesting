from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from pyklip.fakes import retrieve_planet_flux, retrieve_planet
from pyklip.klip import _rotate_wcs_hdr
from copy import copy
from parameter_test_infrastructure import make_dn_per_contrast, FWHMIOWA_calculator, pasep_to_xy
from pyklip.instruments.CHARIS import CHARISData
import pandas as pd
from glob import glob
from astropy.wcs.wcs import FITSFixedWarning
import warnings
import os

# Setting Up (Replicable) Quasi-Random Selection of Files
np.random.seed(314)

warnings.simplefilter('ignore', FITSFixedWarning)
if not os.path.exists('retrievals_companion'):
    os.mkdir('retrievals_companion')


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


orig_filepaths = 'HD1160_cubes/*.fits'
dataset = make_dn_per_contrast(CHARISData(glob(orig_filepaths)))
dn_per_contrast = dataset.dn_per_contrast
rot_angles = dataset.PAs
flipx = dataset.flipx

filepaths = glob('HD1160/klipped_cubes_Wfakes/*')
wavelength_index = 10

files_to_use = [filepaths[np.random.randint(len(filepaths))] for _ in range(20)]
with fits.open(files_to_use[5]) as hdulist:
    output_wcs = WCS(hdulist[0].header, naxis=[1, 2])
    _rotate_wcs_hdr(output_wcs, rot_angles[wavelength_index], flipx=flipx)

df = pd.DataFrame(columns=['Flux', 'FWHM', 'X', 'Y', 'Params'])

for filepath in files_to_use:
    with fits.open(filepath) as hdulist:
        cube = copy(hdulist[1].data)
        dataset_center = [hdulist[1].header['PSFCENTX'], hdulist[1].header['PSFCENTY']]
        dataset_fwhm, _, _ = FWHMIOWA_calculator(hdulist)

    frame = cube[wavelength_index] / dn_per_contrast[wavelength_index]

    sep = 48.33
    pa = 243.93

    fake_flux, xfit, yfit, ret_fwhm = retrieve_planet(frame, dataset_center, output_wcs, sep, pa,
                                                      guessfwhm=dataset_fwhm, guesspeak=1e-4, refinefit=False,
                                                      searchrad=1)

    df = df.append({'Flux': fake_flux, 'FWHM': ret_fwhm, 'X': xfit, 'Y': yfit, 'Params': valuefinder(filepath, 'all')},
                   ignore_index=True)

df.set_index('Params', inplace=True)
df.to_excel(f'retrievals_companion/retrieve_planet_data.xlsx')
