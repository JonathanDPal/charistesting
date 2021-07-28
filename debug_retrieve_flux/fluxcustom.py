if __name__ == '__main__':
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


def retrieve_planet_flux(frame, pa, sep, rotated_wcs, dataset_center, dataset_fwhm, guess_peak_flux=None,
                         force_fwhm=False, searchradius=None, return_all=False):
    def guassian(xy, peak, Fwhm, offset):
        x, y = xy
        sigma = Fwhm / (2 * np.sqrt(2 * np.log(2)))
        return peak * np.exp(-(x ** 2 + y ** 2) / (2 * sigma ** 2)) + offset

    def guassian_force_fwhm(xy_fwhm, peak, offset):
        x, y, Fwhm = xy_fwhm
        sigma = Fwhm / (2 * np.sqrt(2 * np.log(2)))
        return peak * np.exp(-(x ** 2 + y ** 2) / (2 * sigma ** 2)) + offset

    theta = convert_pa_to_image_polar(pa, rotated_wcs)
    x0 = int(round(sep * np.cos(np.radians(theta)) + dataset_center[0]))
    y0 = int(round(sep * np.sin(np.radians(theta)) + dataset_center[1]))

    if searchradius is None:
        searchradius = int(np.ceil(dataset_fwhm))

    preliminary_searchbox = np.copy(frame[y0 - searchradius: y0 + searchradius + 1,
                                          x0 - searchradius: x0 + searchradius + 1])
    preliminary_searchbox[np.where(np.isnan(preliminary_searchbox))] = 0

    # since sometimes the theta found from the position angle and wcs information can be slighly off, we're going to
    # just make sure that we center the Gaussian at the brightest pixel. Usually this will leave x0 and y0
    # unchanged since (x0, y0) will already be the location of the brightest pixel.
    maxindices = np.unravel_index(preliminary_searchbox.argmax(), preliminary_searchbox.shape)
    y0, x0 = maxindices[1] + y0 - searchradius, maxindices[0] + x0 - searchradius

    # now this will be the actual searchbox
    searchbox = np.copy(frame[y0 - searchradius: y0 + searchradius + 1, x0 - searchradius: x0 + searchradius + 1])
    searchbox[np.where(np.isnan(searchbox))] = 0
    uppery, upperx = searchbox.shape

    # building out arrays so that we get a full coordinate system when they are zipped together
    yvals, xvals = np.arange(0, uppery, 1.0) - searchradius, np.arange(0, upperx, 1.0) - searchradius
    numyvals, numxvals = len(yvals), len(xvals)
    yvals = np.array(list(yvals) * numxvals)
    xvals = np.array([[xval] * numyvals for xval in xvals]).flatten()

    vals = [searchbox[int(y + searchradius), int(x + searchradius)] for y, x in zip(yvals, xvals)]

    vals_w_distances = pd.DataFrame({'vals': vals, 'distance squred': [y ** 2 + x ** 2 for y, x in zip(yvals, xvals)]})
    data_to_fit = vals_w_distances[vals_w_distances['distance squared'] <= searchradius ** 2]['vals']

    if guess_peak_flux is None:
        # if None, starting guess (default 1) will be outside of the boundaries we specify, yielding an error
        guess_peak_flux = np.max(data_to_fit) * 0.9

    if force_fwhm:
        guesses = [guess_peak_flux]
        # peak flux cannot be less than zero or greater than the brightest pixel in the image; no bound on offset
        bounds = ((0, -np.inf), (np.max(data_to_fit), np.inf))
    else:
        guesses = [guess_peak_flux, dataset_fwhm]
        bounds = ((0, 0, -np.inf), (np.max(data_to_fit), np.inf, np.inf))  # FWHM cannot be less than zero

    if force_fwhm:
        fwhmlist = np.array([dataset_fwhm] * numxvals * numyvals)
        coordinates = (xvals, yvals, fwhmlist)
        optimalparams, covariance_matrix = curve_fit(f=guassian_force_fwhm, xdata=coordinates, ydata=data_to_fit,
                                                     p0=guesses, bounds=bounds)
    else:
        coordinates = (xvals, yvals)
        optimalparams, covariance_matrix = curve_fit(f=guassian, xdata=coordinates, ydata=data_to_fit, p0=guesses,
                                                     bounds=bounds)

    if return_all:
        return optimalparams
    else:
        return optimalparams[0]  # just the peak flux


if __name__ == '__main__':
    orig_filepaths = 'HD1160_cubes/*.fits'
    dataset = make_dn_per_contrast(CHARISData(glob(orig_filepaths)))
    dn_per_contrast = dataset.dn_per_contrast
    rot_angles = dataset.PAs
    flipx = dataset.flipx

    filepaths = glob('HD1160/klipped_cubes_Wfakes/*')
    filepaths = [filepaths[index] for index in [np.random.randint(len(filepaths)) for _ in range(10)]]

    wavelength_index = 10

    # with fits.open('HD1160/klipped_cubes_Wfakes/HD1160_withfakes_6Annuli_4Subsections_1.0Movement_NoneSpectrum_0'
    #                '.0Smooth_5.0Highpass_-KL50-speccube.fits') as hdulist:
    #

    for filepath in filepaths:
        with fits.open(filepath) as hdulist:
            cube = copy(hdulist[1].data)
            dataset_center = [hdulist[1].header['PSFCENTX'], hdulist[1].header['PSFCENTY']]
            dataset_fwhm, _, _ = FWHMIOWA_calculator(hdulist)
            output_wcs = WCS(hdulist[0].header, naxis=[1, 2])
            _rotate_wcs_hdr(output_wcs, rot_angles[wavelength_index], flipx=flipx)

        frame = cube[wavelength_index] / dn_per_contrast[wavelength_index]

        mask_xy = [144, 80]

        x_pos, y_pos = mask_xy
        ydat, xdat = np.indices(frame.shape)
        distance_from_planet = np.sqrt((xdat - x_pos) ** 2 + (ydat - y_pos) ** 2)
        frame[np.where(distance_from_planet <= 2 * dataset_fwhm)] = np.nan

        fake_fluxes = [5e-4, 5e-5, 5e-6, 1e-4, 1e-5, 1e-6]  # List of Float(s)
        fake_seps = [20, 40, 60]  # List of Integer(s) and/or Float(s)
        fake_PAs = [19, 79, 139, 199, 259, 319]  # List of Integer(s) and/or Float(s)

        retrieved_fluxes, ret_fwhms, ret_xfits, ret_yfits = list(), list(), list(), list()

        for sep in fake_seps:
            for pa in fake_PAs:
                # fake_flux, xfit, yfit, ret_fwhm = retrieve_planet(frame, dataset_center, output_wcs, sep, pa,
                #                                                    guessfwhm=dataset_fwhm, guesspeak=1e-4,
                #                                                    refinefit=False,
                #                                                    searchrad=1)
                # fake_flux = retrieve_planet_flux(frame, pa, sep, output_wcs, dataset_center, dataset_fwhm,
                #                                  force_fwhm=False)
                fake_flux,
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

        # locs = pasep_to_xy(fake_PAs, fake_seps)

        df = pd.DataFrame(index=([str(num+1) for num in range(6)] + ['MEAN']))

        df['Injected Flux (20 sep)'] = at20
        df['Retrieved Flux (20 sep)'] = retrieved_fluxes[:6] + [np.mean(retrieved_fluxes[:6])]
        # df['Injected X (20 sep)'] = [loc[0] for loc in locs[:6]] + [np.nan]
        # df['Retrieved X (20 sep)'] = ret_xfits[:6] + [np.nan]
        # df['Injected Y (20 sep)'] = [loc[1] for loc in locs[:6]] + [np.nan]
        # df['Retrieved Y (20 sep)'] = ret_yfits[:6] + [np.nan]
        # df['Retrieved FWHM (20 sep)'] = ret_fwhms[:6] + [np.mean(ret_fwhms[:6])]

        df['Injected Flux (40 sep)'] = at40
        df['Retrieved Flux (40 sep)'] = retrieved_fluxes[6:12] + [np.mean(retrieved_fluxes[6:12])]
        # df['Injected X (40 sep)'] = [loc[0] for loc in locs[6:12]] + [np.nan]
        # df['Retrieved X (40 sep)'] = ret_xfits[6:12] + [np.nan]
        # df['Injected Y (40 sep)'] = [loc[1] for loc in locs[6:12]] + [np.nan]
        # df['Retrieved Y (40 sep)'] = ret_yfits[6:12] + [np.nan]
        # df['Retrieved FWHM (40 sep)'] = ret_fwhms[6:12] + [np.mean(ret_fwhms[6:12])]

        df['Injected Flux (60 sep)'] = at60
        df['Retrieved Flux (60 sep)'] = retrieved_fluxes[12:] + [np.mean(retrieved_fluxes[12:])]
        # df['Injected X (60 sep)'] = [loc[0] for loc in locs[12:]] + [np.nan]
        # df['Retrieved X (60 sep)'] = ret_xfits[12:] + [np.nan]
        # df['Injected Y (60 sep)'] = [loc[1] for loc in locs[12:]] + [np.nan]
        # df['Retrieved Y (60 sep)'] = ret_yfits[12:] + [np.nan]
        # df['Retrieved FWHM (60 sep)'] = ret_fwhms[12:] + [np.mean(ret_fwhms[12:])]

        df.to_excel(f'retrievals/{valuefinder(filename=filepath, param="all")}retrieve_planet_data.xlsx')
