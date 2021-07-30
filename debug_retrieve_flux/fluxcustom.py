import matplotlib.pyplot as plt
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


def retrieve_planet_flux(frame, pa, sep, output_wcs, dataset_center, dataset_fwhm, theta=None, guess_peak_flux=None,
                         force_fwhm=False, searchradius=None, use_absvals=False, return_all=False, return_r2=False,
                         return_plotting_info=False):
    """
    Identifies the peak flux of a planet by fitting a 2-dimensional Gaussian function.
    ---
    Args:
        frame (2D ndarray): An image from a KLIP output which has been calibrated
        pa (float/int): The position angle of the planet whose flux we are retrieving. Used to calculate theta. If
                        theta is specified, this value will not be used (so you can just input None).
        sep (float/int): The seperation of the planet whose flux we are retreiving.
        output_wcs (astropy.wcs.WCS object): Only output_wcs.cd is what will be used, so this attribute must be
                                             present (and correct).
        dataset_center (list): List of form [x-pos, ypos] specifying the center of the star's PSF.
        dataset_fwhm: The FWHM associated with the observations.
        theta (float): Default: None. The counterclockwise angle in degrees from the x-axis at which the planet is
                       located. If None, will be calculated using pa and output_wcs.
        guess_peak_flux: Default: None. Guess for the peak flux of the planet. If None, will calculate a guess based on
                         the values in the area of the image around the planet.
        force_fwhm (bool): Default: False. If set to True, the curve fitter will be forced to use dataset_fwhm to
                           calculate sigma on the gaussian model. Setting to True is HIGHLY DISCOURAGED because in
                           most cases so far, it has been impossible to fit a gaussian model to the data with this
                           as the FWHM, so usually the scipy.optimize.curve_fit function just silently returns
                           whatever was put in as the guess for all params.
        searchradius (int): Default: None. If None, will use the FWHM as the radius.
        return_all (bool): Default: False. If True, function will return all parameters -- either (peakflux, fwhm,
                           offset) if force_fwhm is False or (peakflux, offset) if force_fwhm is True.
        return_r2 (bool): Default: False. If True, function will return r^2 for the fitted parameters in addition to
                          other returns.
    ---
    Returns:
        Always returns the peak flux. If return_all or return_r2 is True, additional values will be returned.
    """

    def get_r2(actual_vals, predictions):
        ssr = np.sum([(actual_val - prediction) ** 2 for actual_val, prediction in zip(actual_vals, predictions)])
        sst = np.sum([(actual_val - np.mean(actual_vals)) ** 2 for actual_val in actual_vals])
        return 1 - float(ssr) / float(sst)

    def gaussian(xy, peak, Fwhm, offset, y0, x0):
        y, x = xy
        sigma = Fwhm / (2 * np.sqrt(2 * np.log(2)))
        return peak * np.exp(-((y-y0) ** 2 + (x-x0) ** 2) / (2 * sigma ** 2)) + offset

    def gaussian_force_fwhm(xy_fwhm, peak, offset):
        y, x, Fwhm = xy_fwhm
        sigma = Fwhm / (2 * np.sqrt(2 * np.log(2)))
        return peak * np.exp(-(y ** 2 + x ** 2) / (2 * sigma ** 2)) + offset

    if theta is None:
        theta = convert_pa_to_image_polar(pa, output_wcs)

    if searchradius is None:
        searchradius = int(np.ceil(dataset_fwhm))

    y0 = int(round(sep * np.sin(np.radians(theta)) + dataset_center[1]))
    x0 = int(round(sep * np.cos(np.radians(theta)) + dataset_center[0]))

    searchbox = np.copy(frame[y0 - searchradius: y0 + searchradius + 1,
                                    x0 - searchradius: x0 + searchradius + 1])
    searchbox[np.where(np.isnan(searchbox))] = 0

    uppery, upperx = searchbox.shape

    # building out arrays so that we get a full coordinate system when they are zipped together
    yvals, xvals = np.arange(0, uppery, 1.0) - searchradius, np.arange(0, upperx, 1.0) - searchradius
    numyvals, numxvals = len(yvals), len(xvals)
    yvals = np.array(list(yvals) * numxvals)
    xvals = np.array([[xval] * numyvals for xval in xvals]).flatten()

    vals = [searchbox[int(y + searchradius), int(x + searchradius)] for y, x in zip(yvals, xvals)]

    # narrowing down search to just get stuff within the search radius
    vals_w_distances = pd.DataFrame({'vals': vals, 'distance squared': [y ** 2 + x ** 2 for y, x in zip(yvals, xvals)],
                                     'y': yvals, 'x': xvals})
    within_search_radius = vals_w_distances[vals_w_distances['distance squared'] <= searchradius ** 2]
    data_to_fit, yvals, xvals = within_search_radius['vals'], within_search_radius['y'], within_search_radius['x']

    if use_absvals:
        data_to_fit = np.abs(data_to_fit)

    if guess_peak_flux is None:
        # if None, starting guess (default 1) will be outside of the boundaries we specify, yielding an error
        guess_peak_flux = np.max(data_to_fit) * 0.9  # optimal peak flux is usually pretty close to brightest pixel

    if force_fwhm:
        guesses = [guess_peak_flux, 0]
        # peak flux cannot be less than zero or greater than the brightest pixel in the image. upper bound on flux is
        # relatively strict here (as opposed to using star flux as an upper bound, for example) since negative
        # values in the image due to KLIP oversubtraction can cause Gaussian model to be super sharp and overestimate
        # flux by multiple orders of magnitude.
        bounds = ((0, -np.inf), (np.max(data_to_fit), np.inf))  # offset has no lower or upper bound
    else:
        guesses = [guess_peak_flux, dataset_fwhm, 0, 0, 0]
        # FWHM cannot be less than zero; x and y center cannot be adjusted by more than 3 pixels in any direction
        bounds = ((0, 0, -np.inf, -3, -3), (np.max(data_to_fit), np.inf, np.inf, 3, 3))

    if force_fwhm:
        fwhmlist = np.array([dataset_fwhm] * numxvals * numyvals)
        coordinates = (yvals, xvals, fwhmlist)
        optimalparams, covariance_matrix = curve_fit(f=gaussian_force_fwhm, xdata=coordinates, ydata=data_to_fit,
                                                     p0=guesses, bounds=bounds)
    else:
        coordinates = (yvals, xvals)
        optimalparams, covariance_matrix = curve_fit(f=gaussian, xdata=coordinates, ydata=data_to_fit, p0=guesses,
                                                     bounds=bounds)

    if not (return_all or return_r2):  # most of the time this is going to be end of function
        return optimalparams[0]  # just the peak flux

    # this section is just intended to be available for getting more information if needed to assess model & fit
    else:
        if return_r2:
            if force_fwhm:
                predictions = [gaussian_force_fwhm((y, x, fwhm), optimalparams[0], optimalparams[1])
                               for y, x in zip(yvals, xvals)]
            else:
                predictions = [gaussian((y, x), optimalparams[0], optimalparams[1], optimalparams[2],
                                        optimalparams[3], optimalparams[4])
                               for y, x in zip(yvals, xvals)]
            r2 = get_r2(actual_vals=data_to_fit, predictions=predictions)

        if return_all:
            if return_r2:
                if return_plotting_info:
                    return optimalparams, r2, (yvals, xvals), data_to_fit, (float(optimalparams[0]) / float(np.max(
                        data_to_fit)))
                else:
                    return optimalparams, r2
            else:
                return optimalparams
        else:
            return optimalparams[0], r2


def scatter_plot(xy, vals, params, title, absolute, r2):
    def gaussian(xy, peak, fwhm, offset, y0, x0):
        y, x = xy
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        return peak * np.exp(-(((y - y0) ** 2 + x - x0) ** 2) / (2 * sigma ** 2)) + offset

    distances = [y ** 2 + x ** 2 for y, x in zip(xy[0], xy[1])]
    predictions = [gaussian(d, params[0], params[1], params[2], params[3], params[4])
                   for d in np.linspace(0, np.max(distances), 200)]
    fig = plt.figure()
    plt.plot(np.linspace(0, np.max(distances), 200), predictions)
    plt.scatter(distances, vals)
    plt.xlabel('distance')
    if absolute:
        plt.ylabel('absolute value')
        title = f'{title}_absolute'
    else:
        plt.ylabel('value')
    plt.title(title, y=1.07)
    plt.text(0.9, 0.9, rf'$r^{2} = {round(r2, 2)}$', horizontalalignment='center',
             verticalalignment='center', transform=fig.transFigure)
    plt.text(0.5, 0.9, f'Flux: {f"{params[0]:0.3e}"}, FWHM: {round(params[1], 3)}, '
                        f'Offset: {f"{params[2]:0.3e}"}'.replace("'", ""),
             horizontalalignment='center', verticalalignment='center', transform=fig.transFigure)
    plt.savefig(f'retrievals/{title}_scatter.png')
    plt.close(fig)


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
    r2s = dict()
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

        # retrieved_fluxes, ret_fwhms, ret_xfits, ret_yfits = list(), list(), list(), list()
        ret_fluxes, ret_fwhms, ret_offsets, ret_r2s = list(), list(), list(), list()  # ret = retrieved
        a_ret_fluxes, a_ret_fwhms, a_ret_offsets, a_ret_r2s = list(), list(), list(), list()
        ret_percent, a_ret_percent = list(), list()
        for sep in fake_seps:
            for pa in fake_PAs:
                # fake_flux, xfit, yfit, ret_fwhm = retrieve_planet(frame, dataset_center, output_wcs, sep, pa,
                #                                                    guessfwhm=dataset_fwhm, guesspeak=1e-4,
                #                                                    refinefit=False,
                #                                                    searchrad=1)
                # fake_flux = retrieve_planet_flux(frame, pa, sep, output_wcs, dataset_center, dataset_fwhm,
                #                                  force_fwhm=False)
                # for use_abs in [True, False]:
                for use_abs in [False]:
                    optimalparams, r2, xy, data_to_fit, percent = \
                        retrieve_planet_flux(frame, pa, sep, output_wcs, dataset_center, dataset_fwhm, return_all=True,
                                             return_r2=True, return_plotting_info=True, use_absvals=use_abs)
                    # fake_flux, ret_fwhm, offset, ret_x0, ret_y0 = optimalparams
                    params = f'{valuefinder(filepath, "all")}'.replace(' ', '')
                    title = f'{params}_sep-{sep}_pa-{pa}'
                    r2s[title] = r2
                    # scatter_plot(xy=xy, vals=data_to_fit, params=optimalparams, title=title,
                    #              absolute=use_abs, r2=r2)

                    # if use_abs:
                    #     a_ret_fluxes.append(fake_flux)
                    #     a_ret_fwhms.append(ret_fwhm)
                    #     a_ret_offsets.append(offset)
                    #     a_ret_r2s.append(r2)
                    #     a_ret_percent.append(percent)
                    # else:
                    #     ret_fluxes.append(fake_flux)
                    #     ret_fwhms.append(ret_fwhm)
                    #     ret_offsets.append(offset)
                    #     ret_r2s.append(r2)
                    #     ret_percent.append(percent)
                # ret_xfits.append(xfit - dataset_center[0])
                # ret_yfits.append(yfit - dataset_center[0])

        # at20 = [fake_fluxes[0], fake_fluxes[3]] * 3
        # at20.append(float(np.mean(at20)))
        # at40 = [fake_fluxes[1], fake_fluxes[4]] * 3
        # at40.append(float(np.mean(at40)))
        # at60 = [fake_fluxes[2], fake_fluxes[5]] * 3
        # at60.append(float(np.mean(at60)))
        #
        # assert len(ret_fluxes) == 18
        # assert len(a_ret_fluxes) == 18
        #
        # # locs = pasep_to_xy(fake_PAs, fake_seps)
        #
        # df = pd.DataFrame(index=([str(num + 1) for num in range(6)] + ['MEAN']))
        #
        # df['Injected Flux (20 sep)'] = at20
        # df['Retrieved Flux (20 sep)'] = ret_fluxes[:6] + [np.mean(ret_fluxes[:6])]
        # df['% of max no abs (20 sep)'] = ret_percent[:6] + [np.mean(ret_percent[:6])]
        # df['Retrieved Flux w/ abs (20 sep)'] = a_ret_fluxes[:6] + [np.mean(a_ret_fluxes[:6])]
        # df['% of max w/ abs (20 sep)'] = a_ret_percent[:6] + [np.mean(a_ret_percent[:6])]
        # # df['Injected X (20 sep)'] = [loc[0] for loc in locs[:6]] + [np.nan]
        # # df['Retrieved X (20 sep)'] = ret_xfits[:6] + [np.nan]
        # # df['Injected Y (20 sep)'] = [loc[1] for loc in locs[:6]] + [np.nan]
        # # df['Retrieved Y (20 sep)'] = ret_yfits[:6] + [np.nan]
        # df['Retrieved FWHM (20 sep)'] = ret_fwhms[:6] + [np.mean(ret_fwhms[:6])]
        # df['Retrieved FWHM w/ abs (20 sep)'] = a_ret_fwhms[:6] + [np.mean(a_ret_fwhms[:6])]
        # df['Retrieved Offset (20 sep)'] = ret_offsets[:6] + [np.mean(ret_offsets[:6])]
        # df['Retrieved Offset w/ abs (20 sep)'] = a_ret_offsets[:6] + [np.mean(a_ret_offsets[:6])]
        # df['R^2 (20 sep)'] = ret_r2s[:6] + [np.mean(ret_r2s[:6])]
        # df['R^2 w/ abs (20 sep)'] = a_ret_r2s[:6] + [np.mean(a_ret_r2s[:6])]
        #
        # df['Injected Flux (40 sep)'] = at40
        # df['Retrieved Flux (40 sep)'] = ret_fluxes[6:12] + [np.mean(ret_fluxes[6:12])]
        # df['% of max no abs (40 sep)'] = ret_percent[6:12] + [np.mean(ret_percent[6:12])]
        # df['Retrieved Flux w/ abs (40 sep)'] = a_ret_fluxes[6:12] + [np.mean(a_ret_fluxes[6:12])]
        # df['% of max w/ abs (40 sep)'] = a_ret_percent[6:12] + [np.mean(a_ret_percent[6:12])]
        # # df['Injected X (40 sep)'] = [loc[0] for loc in locs[6:12]] + [np.nan]
        # # df['Retrieved X (40 sep)'] = ret_xfits[6:12] + [np.nan]
        # # df['Injected Y (40 sep)'] = [loc[1] for loc in locs[6:12]] + [np.nan]
        # # df['Retrieved Y (40 sep)'] = ret_yfits[6:12] + [np.nan]
        # df['Retrieved FWHM (40 sep)'] = ret_fwhms[6:12] + [np.mean(ret_fwhms[6:12])]
        # df['Retrieved FWHM w/ abs (40 sep)'] = a_ret_fwhms[6:12] + [np.mean(a_ret_fwhms[6:12])]
        # df['Retrieved Offset (40 sep)'] = ret_offsets[6:12] + [np.mean(ret_offsets[6:12])]
        # df['Retrieved Offset w/ abs (40 sep)'] = a_ret_offsets[6:12] + [np.mean(a_ret_offsets[6:12])]
        # df['R^2 (40 sep)'] = ret_r2s[6:12] + [np.mean(ret_r2s[6:12])]
        # df['R^2 w/ abs (40 sep)'] = a_ret_r2s[6:12] + [np.mean(a_ret_r2s[6:12])]
        #
        # df['Injected Flux (60 sep)'] = at60
        # df['Retrieved Flux (60 sep)'] = ret_fluxes[12:] + [np.mean(ret_fluxes[12:])]
        # df['% of max no abs (60 sep)'] = ret_percent[12:] + [np.mean(ret_percent[12:])]
        # df['Retrieved Flux w/ abs (60 sep)'] = a_ret_fluxes[12:] + [np.mean(a_ret_fluxes[12:])]
        # df['% of max w/ abs (60 sep)'] = a_ret_percent[12:] + [np.mean(a_ret_percent[12:])]
        # # df['Injected X (60 sep)'] = [loc[0] for loc in locs[12:]] + [np.nan]
        # # df['Retrieved X (60 sep)'] = ret_xfits[12:] + [np.nan]
        # # df['Injected Y (60 sep)'] = [loc[1] for loc in locs[12:]] + [np.nan]
        # # df['Retrieved Y (60 sep)'] = ret_yfits[12:] + [np.nan]
        # df['Retrieved FWHM (60 sep)'] = ret_fwhms[12:] + [np.mean(ret_fwhms[12:])]
        # df['Retrieved FWHM w/ abs (60 sep)'] = a_ret_fwhms[12:] + [np.mean(a_ret_fwhms[12:])]
        # df['Retrieved Offset (60 sep)'] = ret_offsets[12:] + [np.mean(ret_offsets[12:])]
        # df['Retrieved Offset w/ abs (60 sep)'] = a_ret_offsets[12:] + [np.mean(a_ret_offsets[12:])]
        # df['R^2 (60 sep)'] = ret_r2s[12:] + [np.mean(ret_r2s[12:])]
        # df['R^2 w/ abs (60 sep)'] = a_ret_r2s[12:] + [np.mean(a_ret_r2s[12:])]

        # df.to_excel(f'retrievals/{valuefinder(filename=filepath, param="all")}retrieve_planet_data.xlsx'.replace(' ',
        #                                                                                                          ''))

    for key in r2s.keys():
        with ('rsquared.txt', 'a') as file:
            file.write(f'\n{key}-{r2s[key]}')
