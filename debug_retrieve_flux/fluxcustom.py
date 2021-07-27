from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
# from pyklip.fakes import retrieve_planet_flux, retrieve_planet
from pyklip.fakes import convert_pa_to_image_polar
from pyklip.klip import _rotate_wcs_hdr
from copy import copy
from parameter_test_infrastructure import make_dn_per_contrast, FWHMIOWA_calculator, pasep_to_xy
from pyklip.instruments.CHARIS import CHARISData
import pandas as pd
from glob import glob
from astropy.wcs.wcs import FITSFixedWarning
import warnings
import os
from scipy.optimize import curve_fit

warnings.simplefilter('ignore', FITSFixedWarning)
if not os.path.exists('retrievals'):
    os.mkdir('retrievals')

np.random.seed(314)


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


def retrieve_planet_flux(data_frame, PA, Sep, rotated_wcs, center, fwhm, guess_peak_flux=None, force_fwhm=True,
                         searchrad=None):
    def guassian(xy, a, sigma):
        x, y = xy
        return np.ravel(a * np.exp(-(x ** 2 + y ** 2) / (2 * sigma ** 2)))

    def guassian_force_fwhm(xy_fwhm, a):
        x, y, Fwhm = xy_fwhm
        sigma = Fwhm / (2 * np.sqrt(2 * np.log(2)))
        return np.ravel(a * np.exp(-(x ** 2 + y ** 2) / (2 * sigma ** 2)))

    theta = convert_pa_to_image_polar(PA, rotated_wcs)
    x0 = int(round(Sep * np.cos(np.radians(theta)) + center[0]))
    y0 = int(round(Sep * np.sin(np.radians(theta)) + center[1]))

    if guess_peak_flux is None:
        # if left as None, the initial guess (1) will be outside of the boundaries we specify later, throwing an error
        guess_peak_flux = 1e-3

    if force_fwhm:
        guesses = [guess_peak_flux]
        bounds = (0, 1)  # peak flux cannot be less than zero or greater than star flux
    else:
        guesses = [guess_peak_flux, fwhm]
        bounds = (0, (1, np.inf))  # sigma cannot be less than zero

    if searchrad is None:
        searchrad = int(np.ceil(fwhm))

    searchbox = np.copy(data_frame[y0 - searchrad: y0 + searchrad + 1, x0 - searchrad: x0 + searchrad + 1])
    searchbox[np.where(np.isnan(searchbox))] = 0
    upperx = searchbox.shape[1]
    uppery = searchbox.shape[0]

    if force_fwhm:
        coordinates = np.meshgrid(np.arange(0, upperx, 1.0) - x0, np.arange(0, uppery, 1.0) - y0, fwhm)
        optimalparams, covariance = curve_fit(f=guassian_force_fwhm, xdata=coordinates, ydata=np.ravel(searchbox),
                                              p0=guesses, bounds=bounds)
    else:
        coordinates = np.meshgrid(np.arange(0, upperx, 1.0) - x0, np.arange(0, uppery, 1.0) - y0)
        optimalparams, covariance = curve_fit(f=guassian, xdata=coordinates, ydata=np.ravel(searchbox), p0=guesses,
                                              bounds=bounds)

    peakflux = optimalparams[0]
    return peakflux


orig_filepaths = 'HD1160_cubes/*.fits'
dataset = make_dn_per_contrast(CHARISData(glob(orig_filepaths)))
dn_per_contrast = dataset.dn_per_contrast
rot_angles = dataset.PAs
flipx = dataset.flipx

filepaths = glob('HD1160/klipped_cubes_Wfakes/*')
filepaths = [filepaths[index] for index in [np.random.randint(len(filepaths)) for _ in range(10)]]

wavelength_index = 10

with fits.open('HD1160/klipped_cubes_Wfakes/HD1160_withfakes_6Annuli_4Subsections_1.0Movement_NoneSpectrum_0.0Smooth_5'
               '.0Highpass_-KL50-speccube.fits') as hdulist:
    output_wcs = WCS(hdulist[0].header, naxis=[1, 2])
    _rotate_wcs_hdr(output_wcs, rot_angles[wavelength_index], flipx=flipx)

for filepath in filepaths:
    with fits.open(filepath) as hdulist:
        cube = copy(hdulist[1].data)
        dataset_center = [hdulist[1].header['PSFCENTX'], hdulist[1].header['PSFCENTY']]
        dataset_fwhm, _, _ = FWHMIOWA_calculator(hdulist)

    frame = cube[wavelength_index] / dn_per_contrast[wavelength_index]

    mask_xy = [144, 80]

    x_pos = mask_xy[0]
    y_pos = mask_xy[1]
    ydat, xdat = np.indices(frame.shape)
    distance_from_planet = np.sqrt((xdat - x_pos) ** 2 + (ydat - y_pos) ** 2)
    frame[np.where(distance_from_planet <= 2 * dataset_fwhm)] = np.nan

    fake_fluxes = [5e-4, 5e-5, 5e-6, 1e-4, 1e-5, 1e-6]  # List of Float(s)
    fake_seps = [20, 40, 60]  # List of Integer(s) and/or Float(s)
    fake_PAs = [19, 79, 139, 199, 259, 319]  # List of Integer(s) and/or Float(s)

    retrieved_fluxes = []
    ret_fwhms = []
    ret_xfits = []
    ret_yfits = []
    for sep in fake_seps:
        for pa in fake_PAs:
            # fake_flux, xfit, yfit, ret_fwhm = retrieve_planet(frame, dataset_center, output_wcs, sep, pa,
            #                                                    guessfwhm=dataset_fwhm, guesspeak=1e-4,
            #                                                    refinefit=False,
            #                                                    searchrad=1)
            fake_flux = retrieve_planet_flux(frame, pa, sep, output_wcs, dataset_center, dataset_fwhm)
            retrieved_fluxes.append(fake_flux)
            # ret_fwhms.append(ret_fwhm)
            # ret_xfits.append(xfit - dataset_center[0])
            # ret_yfits.append(yfit - dataset_center[0])

    at20 = [fake_fluxes[0], fake_fluxes[3]] * 3
    at20.append(float(np.mean(at20)))
    at40 = [fake_fluxes[1], fake_fluxes[4]] * 3
    at40.append(float(np.mean(at40)))
    at60 = [fake_fluxes[2], fake_fluxes[5]] * 3
    at60.append(float(np.mean(at60)))

    assert len(retrieved_fluxes) == 18

    act20 = retrieved_fluxes[:6] + [np.mean(retrieved_fluxes[:6])]
    act40 = retrieved_fluxes[6:12] + [np.mean(retrieved_fluxes[6:12])]
    act60 = retrieved_fluxes[12:] + [np.mean(retrieved_fluxes[12:])]

    # fwhms20 = ret_fwhms[:6] + [np.mean(ret_fwhms[:6])]
    # fwhms40 = ret_fwhms[6:12] + [np.mean(ret_fwhms[6:12])]
    # fwhms60 = ret_fwhms[12:] + [np.mean(ret_fwhms[12:])]
    #
    # xfits20 = ret_xfits[:6] + [np.nan]
    # xfits40 = ret_xfits[6:12] + [np.nan]
    # xfits60 = ret_xfits[12:] + [np.nan]
    #
    # yfits20 = ret_yfits[:6] + [np.nan]
    # yfits40 = ret_yfits[6:12] + [np.nan]
    # yfits60 = ret_yfits[12:] + [np.nan]

    # locs = pasep_to_xy(fake_PAs, fake_seps)
    #
    # xlocs20 = [loc[0] for loc in locs[:6]] + [np.nan]
    # xlocs40 = [loc[0] for loc in locs[6:12]] + [np.nan]
    # xlocs60 = [loc[0] for loc in locs[12:]] + [np.nan]
    #
    # ylocs20 = [loc[1] for loc in locs[:6]] + [np.nan]
    # ylocs40 = [loc[1] for loc in locs[6:12]] + [np.nan]
    # ylocs60 = [loc[1] for loc in locs[12:]] + [np.nan]

    df = pd.DataFrame(index=([str(num+1) for num in range(6)] + ['MEAN']))

    df['Injected Flux (20 sep)'] = at20
    df['Retrieved Flux (20 sep)'] = act20
    # df['Injected X (20 sep)'] = xlocs20
    # df['Retrieved X (20 sep)'] = xfits20
    # df['Injected Y (20 sep)'] = ylocs20
    # df['Retrieved Y (20 sep)'] = yfits20
    # df['Retrieved FWHM (20 sep)'] = fwhms20

    df['Injected Flux (40 sep)'] = at40
    df['Retrieved Flux (40 sep)'] = act40
    # df['Injected X (40 sep)'] = xlocs40
    # df['Retrieved X (40 sep)'] = xfits40
    # df['Injected Y (40 sep)'] = ylocs40
    # df['Retrieved Y (40 sep)'] = yfits40
    # df['Retrieved FWHM (40 sep)'] = fwhms40

    df['Injected Flux (60 sep)'] = at60
    df['Retrieved Flux (60 sep)'] = act60
    # df['Injected X (60 sep)'] = xlocs60
    # df['Retrieved X (60 sep)'] = xfits60
    # df['Injected Y (60 sep)'] = ylocs60
    # df['Retrieved Y (60 sep)'] = yfits60
    # df['Retrieved FWHM (60 sep)'] = fwhms60

    df.to_excel(f'retrievals/{valuefinder(filename=filepath, param="all")}retrieve_planet_data.xlsx')
