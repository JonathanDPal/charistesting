from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from glob import glob
from pyklip.instruments.CHARIS import CHARISData
from copy import copy, deepcopy
from pyklip.parallelized import klip_dataset
from pyklip.klip import meas_contrast, _rotate_wcs_hdr
from pyklip.fakes import inject_planet, convert_pa_to_image_polar
from pyklip.fakes import retrieve_planet_flux as pyklip_retrieve_planet_flux
from pyklip.kpp.utils.mathfunc import gauss2d
from pyklip.kpp.metrics.crossCorr import calculate_cc
from pyklip.kpp.stat.statPerPix_utils import get_image_stat_map_perPixMasking
from pyklip.kpp.detection.detection import point_source_detection
import pandas as pd
import sys
import os
from contextlib import contextmanager
import inspect
from multiprocessing.pool import Pool
from time import time
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score


####################
# HELPER FUNCTIONS #
####################

@contextmanager
def log_file_output(directory, write_type='a'):
    """
    Has outputs written out to a log file in the specified directory instead of printed in terminal.
    """
    with open(f'{directory}/log.txt', write_type) as log_file:
        old_stdout = sys.stdout
        sys.stdout = log_file
        try:
            yield
        finally:
            sys.stdout = old_stdout


def FWHMIOWA_calculator(speccubefile, filtname=None):
    """
    This used to calculate FWHM based on the central wavelength, but we found that the equation for that calculation
    was faulty, so at this point, it just spits back 3.5 as the FWHM. Just kept in this format to minimize the
    number of edits needed and also to permit modifications in the future such that we can have different FWHMs for
    different central wavelengths.
    """
    #####
    # If we were to make it formulaic again, it would look something like this:
    # wavelengths = {'j': 1200e-9, 'h': 1550e-9, 'k': 2346e-9, 'broadband': 1550e-9}
    # if filtname is None:
    #     wavelength = wavelengths[str.lower(speccubefile[1].header['FILTNAME'])]
    # else:
    #     wavelength = wavelengths[str.lower(filtname)]
    # D = 8
    # FWHM = (some constant based on instrument characteristics) * wavelength / D
    #########
    # Another way to approach it would just be something like this (with wavelength coming from above stuff):
    # fwhms = {'j': {fwhm1}, 'h': {fwhm2}, 'k':{fwhm3}, 'broadband': 3.5}
    # FWHM = fwhms[wavelength]
    ########
    lenslet_scale = 0.0162
    field_radius = 1.035
    FWHM = 3.5
    IWA = 5
    OWA = (field_radius / lenslet_scale) - FWHM

    return FWHM, IWA, OWA


def make_dn_per_contrast(dataset):
    """
    Calculates and sets spot_ratio and dn_per_contrast attributes for an initialized CHARISData dataset.

    Returns modified CHARISData dataset object.
    """

    # Gets number of input fits files (Ncubes) and number of wavelengths (Nwln)
    Nframes = dataset.input.shape[0]  # This dimension is Ncubes*Nwln
    Ncubes = np.size(np.unique(dataset.filenums))
    Nwln = int(Nframes / Ncubes)

    # Gets wavelength in microns; 1D array with shape (Nfiles * Nwlns,)
    wln_um = dataset.wvs

    # Calculates the spot/star ratio for each wavelength, in um; 1D array with shape (Nfiles * Nwlns,)
    dataset.spot_ratio = 2.72e-3 * (wln_um / 1.55) ** -2

    # Finds the mean spot flux across all files at each wavelength; 1D array with shape (Nwlns,)
    mean_spot_flux = np.nanmean(dataset.spot_fluxes.reshape(Ncubes, Nwln), axis=0)

    # Tiles the mean spot flux array to repeat Ncubes times; 1D array with shape (Nfiles * Nwlns,)
    mean_spot_fluxes = np.tile(mean_spot_flux, Ncubes)

    # Calculates and sets the dn_per_contrast
    dataset.dn_per_contrast = mean_spot_fluxes / dataset.spot_ratio


def pasep_to_xy(PAs, seps):
    """
    Takes lists of position angles and seperations and yields a numpy array with x-y coordinates for each combo, using
    the convention used in the table of values outputted by the planet detection software. The origin used in their
    convention is the center of the star's PSF.
    """
    PAs = [float(pa) for pa in PAs]  # if not a float, then pa / 180 will yield zero in many cases
    radians = np.array(PAs) / 180 * np.pi
    locs = []
    for sep in seps:
        for rad in radians:
            x = -np.sin(rad) * sep
            y = np.cos(rad) * sep
            loc = [x, y]
            locs.append(loc)
    return locs


def distance(xy1, xy2):
    """
    Inputs should be of form [x-coor, y-coor] (list, numpy array, or tuple)
    """
    return np.sqrt((xy1[0] - xy2[0]) ** 2 + (xy1[1] - xy2[1]) ** 2)


def contrast_measurement(trial_string):
    """
    Would be used for parallelization on contrast measurements. Not currently in use since some jankiness on Condor
    was causing enough errors with contrast that parallelizaion on contrast was leaving behind loose daemons on the
    system. Will probably be reintroduced later on
    """
    t = Trial.from_string(trial_string)
    t.get_contrast()


def planet_detection(trial_string):
    """
    Used for parallelization on planet detections. Not currently in use because the planet detection section of
    pyKLIP utilizes the multiprocessing process pool for computing SNR maps and multiprocessing cannot support
    a pool inside of another pool. (I might utilize a workaround on this later on)
    """
    t = Trial.from_string(trial_string)
    t.detect_planets()


def parameter_set_batcher(batchindex, batchsize, args):
    """
    Designed to take a set of parameters and return a subset of combinations so that parameters can be broken up
    into batches using the command line arguments of the parameterprobing scripts.

    Args:
        batchindex: Which batch should be used (USING 1-BASED INDEXING, NOT 0-BASED INDEXING)
        batchsize: How many trials should be in the batches
        args: The KLIP parameters (annuli, subsections, movement, spectrum, corr_smooth, highpass)

    Returns a list or list of lists of parameters to be passed in to TestDataset.
    """
    num_param_combos = np.prod([len(arg) for arg in args])
    remainder = num_param_combos % batchsize
    partial_batch = False  # default value, might be flipped by next section
    if remainder != 0:
        num_full_size_batches = int(np.floor(num_param_combos / batchsize))
        if batchindex > num_full_size_batches:
            partial_batch = True

    annuli, subsections, movement, spectrum, corr_smooth, highpass = args
    params = []
    for ani in annuli:
        for subsec in subsections:
            for mov in movement:
                for spec in spectrum:
                    for cs in corr_smooth:
                        for hp in highpass:
                            params.append((ani, subsec, mov, spec, cs, hp))

    startingindex = batchsize * (batchindex - 1)  # (batchindex - 1) is because it is 1-Based indexing being passed in
    if not partial_batch:
        finalindex = startingindex + batchsize
    else:
        finalindex = None

    newparams = params[startingindex:finalindex]

    return newparams


def retrieve_planet_flux(frame, pa, sep, output_wcs, dataset_center, dataset_fwhm, theta=None, guess_peak_flux=None,
                         force_fwhm=False, searchradius=None):
    """
    Identifies the peak flux of a planet by fitting a 2-dimensional Gaussian function.
    ---
    Required Args:
        frame (2D ndarray): An image from a KLIP output which has been calibrated
        pa (float/int): The position angle of the planet whose flux we are retrieving. Used to calculate theta. If
                        theta is specified, this value will not be used (so you can just input None).
        sep (float/int): The seperation of the planet whose flux we are retreiving.
        output_wcs (astropy.wcs.WCS object): Only output_wcs.cd is what will be used, so this attribute must be
                                             present (and correct).
        dataset_center (list): List of form [x-pos, ypos] specifying the center of the star's PSF.
        dataset_fwhm (float): The FWHM associated with the observations.
    Optional Args:
        theta (float): Default: None. The counterclockwise angle in degrees from the x-axis at which the planet is
                       located. If None, will be calculated using pa and output_wcs.
        guess_peak_flux (float): Default: None. Guess for the peak flux of the planet. If None, will calculate a guess
                                 based on the values in the area of the image around the planet.
        force_fwhm (bool): Default: False. If set to True, the curve fitter will be forced to use dataset_fwhm to
                           calculate sigma on the gaussian model. Setting to True is HIGHLY DISCOURAGED because in
                           most cases so far, it has been impossible to fit a Gaussian model to the data with this
                           as the FWHM, so the scipy.optimize.curve_fit function will just silently return
                           whatever was put in as the guess for all params.
        searchradius (int): Default: None. If None, will use np.ceil(dataset_fwhm) as the search radius.
    ---
    Returns:
        Returns peak flux. If all values within search radius are negative, then the function will return 0.
    """

    def gaussian(xy, peak, Fwhm, offset, y0, x0):
        y, x = xy
        sigma = Fwhm / (2 * np.sqrt(2 * np.log(2)))
        return peak * np.exp(-((y - y0) ** 2 + (x - x0) ** 2) / (2 * sigma ** 2)) + offset

    def gaussian_force_fwhm(xy_fwhm, peak, offset, y0, x0):
        y, x, Fwhm = xy_fwhm
        sigma = Fwhm / (2 * np.sqrt(2 * np.log(2)))
        return peak * np.exp(-((y - y0) ** 2 + (x - x0) ** 2) / (2 * sigma ** 2)) + offset

    if theta is None:
        theta = convert_pa_to_image_polar(pa, output_wcs)

    if searchradius is None:
        searchradius = int(np.ceil(dataset_fwhm))

    y0 = int(round(sep * np.sin(np.deg2rad(theta)) + dataset_center[1]))
    x0 = int(round(sep * np.cos(np.deg2rad(theta)) + dataset_center[0]))

    searchbox = np.copy(frame[y0 - searchradius: y0 + searchradius + 1, x0 - searchradius: x0 + searchradius + 1])
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

    allneg = (np.max(data_to_fit) <= 0)  # boolean value; using <= instead of < since code will break if both upper
    # bound and lower bound are 0 for flux

    if allneg:
        # just making this as a definition (if all negative, then contrast sucks, so nothing can be detected). We're
        # defining this as zero since we divide by this to get calibrated contrast, so this would be contrast=infinity
        return 0
    else:
        if guess_peak_flux is None:
            # if None, starting guess (default 1) will be outside of the boundaries we specify, yielding an error
            guess_peak_flux = np.max(data_to_fit) * 0.9  # optimal peak flux is usually pretty close to brightest pixel,
            # but we don't want to have our guess be all the way at the upper bound (the brightest pixel), so using 90%

        if force_fwhm:
            guesses = [guess_peak_flux, 0]
            # peak flux cannot be less than zero or greater than the brightest pixel in the image. upper bound on
            # flux is relatively strict here (as opposed to using star flux as an upper bound, for example) since
            # negative values in the image due to KLIP oversubtraction can cause Gaussian model to be super sharp and
            # overestimate flux by multiple orders of magnitude.
            bounds = ((0, -np.inf), (np.max(data_to_fit), np.inf))  # offset has no lower or upper bound
        else:
            guesses = [guess_peak_flux, dataset_fwhm, 0, 0, 0]  # zeros are for offset, y0, and x0
            # FWHM cannot be less than zero; x and y center cannot be adjusted by more than 3 pixels in any direction
            bounds = ((0, 0, -np.inf, -3, -3), (np.max(data_to_fit), np.inf, np.inf, 3, 3))

        if force_fwhm:
            fwhmlist = np.array([dataset_fwhm] * numxvals * numyvals)
            coordinates = (yvals, xvals, fwhmlist)
            optimalparams, covariance_matrix = curve_fit(f=gaussian_force_fwhm, xdata=coordinates, ydata=data_to_fit,
                                                         p0=guesses, bounds=bounds)
            # new stuff for testing
            predictions = [gaussian_force_fwhm((y, x, fwhm), optimalparams[0], optimalparams[1], optimalparams[2],
                                               optimalparams[3]) for y, x, fwhm in
                                               zip(coordinates[0], coordinates[1], coordinates[2])]
            r2 = r2_score(data_to_fit, predictions)
        else:
            coordinates = (yvals, xvals)
            optimalparams, covariance_matrix = curve_fit(f=gaussian, xdata=coordinates, ydata=data_to_fit, p0=guesses,
                                                         bounds=bounds)
            # new stuff for testing
            predictions = [gaussian((y, x), optimalparams[0], optimalparams[1], optimalparams[2],
                                            optimalparams[3], optimalparams[4]) for y, x in
                                            zip(coordinates[0], coordinates[1])]
            r2 = r2_score(data_to_fit, predictions)

        return optimalparams, r2  # peak flux


def params_from_text_file(paramsfile):
    """
    Designed to be utilized on a text file generated by check.py (inside organizational_scripts).
    """
    [annuli, subsections, movement, spectrum, numbasis, corr_smooth, highpass] = [list() for _ in range(7)]
    with open(paramsfile) as file:
        lines = [line[:-1] for line in file]  # last character is '\n', which we don't want
    for line in lines:
        [ann, sbs, mov, spec, nb, cs, hp] = [p.replace(' ', '') for p in line.split(',')]
        ann, sbs, nb = int(ann), int(sbs), int(nb)
        mov, cs = float(mov), float(cs)
        if spec.lower() == 'none':
            spec = None
        if hp.lower() == 'true':
            hp = True
        elif hp.lower() == 'false':
            hp = False
        else:
            hp = float(hp)
        annuli.append(ann)
        subsections.append(sbs)
        movement.append(mov)
        spectrum.append(spec)
        numbasis.append(nb)
        corr_smooth.append(cs)
        highpass.append(hp)

    return annuli, subsections, movement, spectrum, numbasis, corr_smooth, highpass


def append_dataset_info(info_filename, dataset_placeholder):
    """
    Adds in information to TestDataset object needed for contrast and detection from a text file. Avoids the need to
    build CHARISData object (memory-intensive). File with information should be an output file from
    'organizational_scripts/get_dataset_info.py'.
    ---
    Params:
        info_filename (str): Name of the file which has the needed information from the CHARISData object
        testdataset (TestDataset): TestDataset object to have information appended to.
    """
    with open(info_filename) as f:
        lines = [line for line in f]
    angles, dn_per_contrast, wvs = list(), list(), list()
    for index, line in enumerate(lines):
        if 'Length' in line:
            dataset_placeholder.leNgth = int(line.split(' ')[1])
            lengthline = index
        elif 'Flip_x' in line:
            dataset_placeholder.flipx = bool(line.split(' ')[1])
            flipxline = index
        elif 'DN_per_contrast' in line:
            dn_per_contrast_start = index + 1
        elif 'Wavelengths' in line:
            wvs_start = index + 1
    linenumber = 1
    while linenumber < flipxline:
        angles.append(float(lines[linenumber]))
        linenumber += 1
    linenumber = dn_per_contrast_start
    while linenumber < wvs_start - 1:
        dn_per_contrast.append(float(lines[linenumber]))
        linenumber += 1
    linenumber = wvs_start
    while linenumber < lengthline:
        wvs.append(float(lines[linenumber]))
    dataset_placeholder.angles = np.array(angles)
    dataset_placeholder.dn_per_contrast = np.array(dn_per_contrast)
    dataset_placeholder.wvs = np.array(wvs)


####################################################################################
# TestDataset Will Manage a List of Trials (one for each group of KLIP Parameters) #
####################################################################################
class Trial:
    """
    NOTE: The user will almost certainly not interact with this class directly, rather they will interact with an
    instance of TestDataset and that instance of TestDataset will interact with instances of this class.
    ---
    Stores a particular set of KLIP parameters and then is able to run contrast measurement or planet detection code
    from KLIP for the KLIP output with that particular set of parameters.
    """

    def __init__(self, object_name, mask_xy, annuli, subsections, movement, numbasis, spectrum, corr_smooth,
                 fake_PAs, fake_fluxes, fake_fwhm, fake_seps, rot_angs, flipx, dn_per_contrast, wln_um, highpass,
                 length):
        self.object_name = object_name
        self.mask_xy = mask_xy

        # for when rebuilding from string
        movement = float(movement)
        corr_smooth = float(corr_smooth)

        # KLIP params used (all are specific values except for numbasis which is still a list of values at this point)
        # By specific, I mean, for example, self.annuli=9 as opposed to self.annuli=[3, 5, 7, 9, 11], like you would
        # see in the self.annuli attribute of TestDataset
        self.annuli = annuli
        self.subsections = subsections
        self.movement = movement
        self.numbasis = numbasis
        self.spectrum = spectrum
        self.corr_smooth = corr_smooth

        self.fake_PAs = fake_PAs
        self.fake_fluxes = fake_fluxes
        self.fake_fwhm = fake_fwhm
        self.fake_seps = fake_seps

        # Needed to Replicate Rotation for Fake Planet Retrieval
        self.rot_angs = rot_angs
        self.flipx = flipx

        self.dn_per_cotrast = np.array(dn_per_contrast)
        self.wln_um = wln_um  # only being used to identify wavelength in output filepath name

        # Switching Highpass To Image Space If Necessary
        if not isinstance(highpass, bool):
            highpass = float(highpass)
            self.highpass = length / (highpass * 2 * np.sqrt(2 * np.log(2)))
        else:
            self.highpass = highpass

        # String Identifying Parameters Used (Used Later For Saving Contrast Info)
        self.klip_parameters = str(annuli) + 'Annuli_' + str(subsections) + 'Subsections_' + str(movement) + \
                               'Movement_' + str(spectrum) + 'Spectrum_' + str(corr_smooth) + 'Smooth_' + str(
            highpass) + 'Highpass_'

        # Filepaths to KLIPped Datacubes
        self.filepaths_Wfakes = [self.object_name + '/klipped_cubes_Wfakes/' + self.object_name + '_withfakes_' +
                                 self.klip_parameters + f'-KL{nb}-speccube.fits' for nb in self.numbasis]
        self.filepaths_Nfakes = [self.object_name + '/klipped_cubes_Nfakes/' + self.object_name + '_withoutfakes_' +
                                 self.klip_parameters + f'-KL{nb}-speccube.fits' for nb in self.numbasis]

        # Filepath to Save Planet Detection Output To
        self.filepath_detections_prefixes = [self.object_name + f'/detections/{self.klip_parameters}_KL{nb}_SNR-'
                                             for nb in self.numbasis]

        # Building a String Which Contains All of the Information For Rebuilding the Trial Instance (this gets used for
        # parallelization)
        params = [object_name, mask_xy, annuli, subsections, movement, numbasis, spectrum, corr_smooth, fake_PAs,
                  fake_fluxes, fake_fwhm, fake_seps, rot_angs, flipx, dn_per_contrast, wln_um, highpass, length]
        modifiedparams = []
        for i in range(len(params)):
            array = False
            if type(params[i]) == np.ndarray:
                if np.ndim(params[i]) == 2:
                    params[i] = [list(subarray) for subarray in params[i]]
                elif np.ndim(params[i]) == 1:
                    params[i] = list(params[i])
                else:
                    raise Warning(f'If you are going to run parallelized contrast, then it will fail due to having a '
                                  f'3+ dimensional numpy array as one of the arguments. The problematic argument is '
                                  f'{params[i]}')
                array = True
            if type(params[i]) == list:
                list_in_list = []
                for j in range(len(params[i])):
                    if type(params[i][j]) == list:
                        m = '[!'
                        for p in params[i][j]:
                            m += f'{p}!'
                        m += ']'
                        list_in_list.append(m)
                    else:
                        list_in_list.append(params[i][j])
                modifiedparams.append(list_in_list)
            else:
                modifiedparams.append(params[i])
            if array:
                params[i] = np.array(params[i])
        self.rebuild_string = '|'.join([str(modifiedparam) for modifiedparam in modifiedparams])

    @staticmethod
    def list_rebuilder(s):
        """
        Takes the rebuild string created by the __init__ method of Trial and returns the original set of parameters
        which was used to build it. Intended to support the class method from_string.
        """
        s = s.replace(' ', '')
        original_params = s.split('|')
        for i in range(len(original_params)):
            pt = original_params[i]
            if pt[0] == '[':
                sub_pts = pt.split(',')
                sub_pts[0] = sub_pts[0][1:]
                sub_pts[-1] = sub_pts[-1][:-1]
                for j in range(len(sub_pts)):
                    sub_pt = sub_pts[j]
                    if sub_pt[0] == "'":
                        sub_sub_pts = sub_pt.split('!')
                        sub_sub_pts = sub_sub_pts[1:-1]
                        for k in range(len(sub_sub_pts)):
                            sspt = sub_sub_pts[k]
                            try:
                                sspt = float(sspt)
                                if int(sspt) == sspt:
                                    sspt = int(sspt)
                            except ValueError:
                                if str.lower(sspt) == 'none':
                                    sspt = None
                                elif str.lower(sspt) == 'true':
                                    sspt = True
                                elif str.lower(sspt) == 'false':
                                    sspt = False
                            finally:
                                sub_sub_pts[k] = sspt
                        sub_pts[j] = sub_sub_pts
                    else:
                        try:
                            sub_pt = float(sub_pt)
                            if int(sub_pt) == sub_pt:
                                sub_pt = int(sub_pt)
                        except ValueError:
                            if str.lower(sub_pt) == 'none':
                                sub_pt = None
                            elif str.lower(sub_pt) == 'true':
                                sub_pt = True
                            elif str.lower(sub_pt) == 'false':
                                sub_pt = False
                        finally:
                            sub_pts[j] = sub_pt
                    original_params[i] = sub_pts
            else:
                try:
                    pt = float(pt)
                    if int(pt) == pt:
                        pt = int(pt)
                except ValueError:
                    if str.lower(pt) == 'none':
                        pt = None
                    elif str.lower(pt) == 'true':
                        pt = True
                    elif str.lower(pt) == 'false':
                        pt = False
                finally:
                    original_params[i] = pt

        return original_params

    @classmethod
    def from_string(cls, rebuild_string):
        """
        Uses the Trial rebuild string from the __init__ method to recreate the Trial object. This is intended so that
        contrast and planet detection can be parallelized.
        """
        p = cls.list_rebuilder(rebuild_string)
        if len(p) != 18:
            raise ValueError("Incorrect number of arguments given for building class Trial using from_string method.")
        return cls(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14],
                   p[15], p[16], p[17])

    def get_contrast(self, contains_fakes=True, overwrite=False):
        """
        Measures contrast at a particular wavelength, then saves contrast data as a CSV. CSV file is saved from a
        Pandas DataFrame, which makes it easy to load the data back into a Pandas DataFrame later on for analysis.
        ---
        Args:
            contains_fakes (bool): Default: True. Whether to use data with fakes (True) or without fakes (False).
            overwrite (bool): Default: False. Whether or not to override filepath if output filepath already exists.
        """
        if contains_fakes:
            filepaths = self.filepaths_Wfakes
        else:
            filepaths = self.filepaths_Nfakes

        # Measuring Contrast For Each Set of KL Modes
        for filepath_index, filepath in enumerate(filepaths):  # filepath_index used to identify number of KL modes
            try:
                with fits.open(filepath) as hdulist:
                    cube = copy(hdulist[1].data)
                    dataset_center = [hdulist[1].header['PSFCENTX'], hdulist[1].header['PSFCENTY']]
                    dataset_fwhm, dataset_iwa, dataset_owa = FWHMIOWA_calculator(hdulist)
                    output_wcs = WCS(hdulist[0].header, naxis=[1, 2])
                corrupt_file = False
            except OSError:
                corrupt_file = True

            if corrupt_file:
                with open(f'{self.object_name}/corrupt_fits_files.txt', 'a') as f:
                    f.write(f'{filepath}\n')
                continue

            for wavelength_index in range(cube.shape[0]):  # making measurements at every wavelength
                # if not overwriting, then check to see if file exists
                wavelength = round(self.wln_um[wavelength_index], 2)  # in microns
                uncal_contrast_output_filepath = self.object_name + f'/uncalibrated_contrast' \
                                                                    f'/{self.klip_parameters}_KL' \
                                                                    f'{self.numbasis[filepath_index]}' \
                                                                    f'_{wavelength}um_contrast.csv'
                cal_contrast_output_filepath = self.object_name + f'/calibrated_contrast/{self.klip_parameters}_KL' \
                                                                  f'{self.numbasis[filepath_index]}' \
                                                                  f'_{wavelength}um_contrast.csv'

                if not overwrite:
                    if os.path.exists(uncal_contrast_output_filepath) and os.path.exists(cal_contrast_output_filepath):
                        continue

                # need to rotate WCS so that we are looking in right spot using the pyKLIP function for it
                local_output_wcs = deepcopy(output_wcs)
                _rotate_wcs_hdr(local_output_wcs, self.rot_angs[wavelength_index], flipx=self.flipx)

                # Taking Slice of Cube and Calibrating It
                frame = cube[wavelength_index] / self.dn_per_contrast[wavelength_index]

                # Applying Mask to Science Target If Location Specified
                if isinstance(self.mask_xy, (list, tuple)):
                    if not isinstance(self.mask_xy[0], (list, tuple)):
                        x_pos = self.mask_xy[0]
                        y_pos = self.mask_xy[1]

                        ydat, xdat = np.indices(frame.shape)
                        distance_from_planet = np.sqrt((xdat - x_pos) ** 2 + (ydat - y_pos) ** 2)
                        frame[np.where(distance_from_planet <= 2 * dataset_fwhm)] = np.nan
                    else:
                        for position in self.mask_xy:
                            x_pos = position[0]
                            y_pos = position[1]

                            ydat, xdat = np.indices(frame.shape)
                            distance_from_planet = np.sqrt((xdat - x_pos) ** 2 + (ydat - y_pos) ** 2)
                            frame[np.where(distance_from_planet <= 2 * dataset_fwhm)] = np.nan

                # Measuring Algorithm Throughput BEFORE Fake Planets Get Masked Out
                if contains_fakes:
                    retrieved_fluxes = []
                    if np.ndim(self.fake_PAs) == 1:
                        pas = self.fake_PAs
                        fluxes = self.fake_fluxes
                    else:  # i.e., multiple tiers
                        pas = self.fake_PAs[0]  # using highest tier planets for calibration)
                        fluxes = self.fake_fluxes[0]
                    optimalparams = list()
                    r2s = list()
                    for sep in self.fake_seps:
                        for pa in pas:
                            try:
                                opparams, r2 = retrieve_planet_flux(frame, pa, sep, local_output_wcs, dataset_center,
                                                                 dataset_fwhm)
                            except RuntimeError:  # if scipy can't find a good model
                                opparams, r2 = [None, None, None, None, None], None
                            optimalparams.append(list(opparams))
                            r2s.append(r2)
                spl = str.split(filepath, '-')
                for s in spl:
                    if 'KL' in s:
                        nb = int(s.replace('KL', ''))
                        break
                nbs = [nb] * len(r2s)
                if len(optimalparams[0]) == 5:
                    columns = 'peak, fwhm, offset, x0, y0, nb, r2\n'
                elif len(optimalparams[0]) == 4:
                    columns = 'peak, offset, x0, y0, nb, r2\n'
                #TEMPORARY#
                print(optimalparams)
                print(self.klip_parameters)
                #TEMPORARY#
                if not os.path.exists(self.object_name + '/fitmeas'):
                    os.mkdir(self.object_name + '/fitmeas')
                with open(self.object_name + f'/fitmeas/fluxretrivaltesting-{self.klip_parameters}.csv', 'w') as f:
                    f.write(columns)
                    for op, r2, nb in zip(optimalparams, r2s, nbs):
                        opstr = f'{op}'.replace('[', '').replace(']', '').replace(' ', '')
                        f.write(opstr + f'{nb},{r2}\n')

    def detect_planets(self, SNR_threshold=2, datasetwithfakes=True, override=False):
        """
        Looks at a KLIPped dataset with fakes and indicates potential planets. Identifies
        ---
        Args:
            SNR_threshold: Default: 2. Set this to the lowest value to be looked at.
            datasetwithfakes (Bool): Default: True. If True, then run planet detection on dataset containing
                                     injected planets; if False, then run planet detection on dataset not containing
                                     injected planets.
            override (bool): Default: False. Whether or not to override filepath if output filepath already exists.
		"""
        if datasetwithfakes:
            filepaths = self.filepaths_Wfakes
        else:
            filepaths = self.filepaths_Nfakes

        for filepath_index, filepath in enumerate(filepaths):  # filepath_index used to identify number of KL modes

            output_filepath = f'{self.filepath_detections_prefixes[filepath_index]}{SNR_threshold}.csv'

            if not override:
                if os.path.exists(output_filepath):
                    continue

            # Actual Start of Process
            try:
                with fits.open(filepath) as hdulist:
                    image = copy(hdulist[1].data)
                    center = [hdulist[1].header['PSFCENTX'], hdulist[1].header['PSFCENTY']]
                corrupt_file = False
            except OSError:  # occurs if the file is corrupt or empty
                corrupt_file = True

            if corrupt_file:
                with open(f'{self.object_name}/corrupt_fits_files.txt', 'a') as f:
                    f.write(f'{filepath}\n')
                continue

            x_grid, y_grid = np.meshgrid(np.arange(-10, 10), np.arange(-10, 10))
            kernel_gauss = gauss2d(x_grid, y_grid)

            # flat spectrum given here for generating cross-coorelated image so that pyKLIP collapses it into one
            # image, instead of giving seperate images for each wavelength
            image_cc = calculate_cc(image, kernel_gauss, spectrum=np.ones(len(image[0])), nans2zero=True)

            SNR_map = get_image_stat_map_perPixMasking(image_cc, centroid=center, mask_radius=5, Dr=2, type='SNR')

            candidates_table = point_source_detection(SNR_map, center, SNR_threshold, pix2as=1, mask_radius=15,
                                                      maskout_edge=10, IWA=None, OWA=None)

            candidates = pd.DataFrame(candidates_table, columns=['Index', 'SNR Value', 'PA', 'Sep (pix)',
                                                                 'Sep (as)', 'x', 'y', 'row', 'col'])

            fakelocs = pasep_to_xy(self.fake_PAs.flatten(), self.fake_seps.flatten())  # where planets were injected

            candidate_locations = zip(candidates['x'], candidates['y'])  # where stuff was detected

            if self.mask_xy is None:
                self.mask_xy = [[1000, 1000]]  # nothing gets identified as science target
            elif not isinstance(self.mask_xy[0], (list, tuple)):
                self.mask_xy = [self.mask_xy]  # making it a list of a list so that it can get iterated over
                # properly

            distances_from_fakes = []  # going to be an additional column of candidates DataFrame
            distances_from_targets = []  # going to be an additional column of candidates DataFrame
            for c in candidate_locations:
                distances = []
                for fl in fakelocs:
                    distances.append(distance(c, fl))
                distances_from_fakes.append(np.min(distances))
                distances2 = []
                for mask in self.mask_xy:
                    mask = np.array(mask) - np.array(center)  # aligning coordinate systems
                    distances2.append(distance(c, mask))
                distances_from_targets.append(np.min(distances2))

            injected = []  # going to be an additional column of candidates DataFrame
            for d1, d2 in zip(distances_from_targets, distances_from_fakes):
                if d1 < self.fake_fwhm * 1.5:
                    injected.append("Science Target")
                elif d2 < self.fake_fwhm * 1.5:
                    injected.append(True)
                else:
                    injected.append(False)

            # appending more information to output for analysis later on
            candidates['Distance From Fakes'] = distances_from_fakes
            candidates['Distance From Targets'] = distances_from_targets
            candidates['Injected'] = injected

            candidates.to_csv(output_filepath, index=False)

    def __eq__(self, other):
        """
        Checks to see if two Trials have the same attributes. Intended for testing out code functionality.
        """
        equal_attributes = list()
        for i, j in zip(inspect.getmembers(self), inspect.getmembers(other)):
            if i[0].startswith('_') or inspect.ismethod(i[1]) or i[0] == 'rebuild_string':
                continue
            else:
                try:
                    equal_attributes.append(i[1] == j[1])
                    if i[1] != j[1]:
                        print(i[0], "\nself: ", i[1], "\nother: ", j[1])
                except ValueError:
                    same = all(i[1]) == all(j[1])
                    equal_attributes.append(same)
                    if not same:
                        print(i[0], "\nself: ", i[1], "\nother: ", j[1])
        for i in range(len(equal_attributes)):
            if isinstance(equal_attributes[i], (list, np.ndarray)):
                equal_attributes[i] = np.sum(equal_attributes[i]) == len(equal_attributes[i])
        return np.sum(equal_attributes) == len(equal_attributes)


#################################################################################################
# Observation Set (eg. HD1160, BetaPic) Will Have An Instance of TestDataset Associated With It #
#################################################################################################
class Dataset_PlaceHolder:
    """
    Helper class for TestDataset
    """
    def __init__(self):
        pass


class TestDataset:
    """
    The main object which the user will interact with. Will load in CHARIS fileset into CHARISData class (see
    pyklip.instruments.CHARIS) and then create an instance of Trial for each set of KLIP parameters to be looked at.
    """

    def __init__(self, fileset, object_name, mask_xy, fake_fluxes, fake_seps, annuli, subsections, movement,
                 numbasis, corr_smooth, highpass, spectrum, fake_fwhm, fake_PAs, mode, batched, overwrite,
                 memorylite, build_all_combos, build_charis_data):
        self.object_name = object_name
        self.mask_xy = mask_xy

        if not os.path.exists(self.object_name):
            try:
                os.mkdir(self.object_name)
            except FileExistsError:
                # unnecessary 99% of the time, but once in a blue moon, I run into issues with this on Condor
                pass

        self.write_to_log(f'Title for Set: {object_name}', 'w')
        self.write_to_log(f'\nFileset: {fileset}')

        if build_all_combos:
            param_names = ['Annuli', 'Subsections', 'Movement', 'Numbasis', 'Corr_Smooth', 'Highpass', 'Spectrum',
                           'Mode', 'Fake Fluxes', 'Fake Seps', 'Fake PAs', 'Fake FWHM']
            params = [annuli, subsections, movement, numbasis, corr_smooth, highpass, spectrum]
            number_of_paramcombos = np.prod([len(p) for p in params])
            self.write_to_log(f'\nNumber of Parameter Combinations: {number_of_paramcombos}')

            for param in [mode, fake_fluxes, fake_seps, fake_PAs, fake_fwhm]:
                params.append(param)
            for name, param in zip(param_names, params):
                self.write_to_log(f'\n{name}: {param}')
        else:  # don't want to write in all values if using text file (could be tens of thousands)
            param_names = ['KLIP Parameters:' 'Mode', 'Fake Fluxes', 'Fake Seps', 'Fake PAs', 'Fake FWHM']
            params = ['from a text file', mode, fake_fluxes, fake_seps, fake_PAs, fake_fwhm]
            for name, param in zip(param_names, params):
                self.write_to_log(f'\n{name}: {param}')

        self.write_to_log_and_print(f'############### STARTING WORK ON {self.object_name} ################\n')

        if build_charis_data == 'true' or build_charis_data == 'temporary':
            with log_file_output(self.object_name):
                self.dataset = CHARISData(glob(fileset))
                make_dn_per_contrast(self.dataset)
                self.dataset.leNgth = self.dataset.input.shape[1]
            self.write_to_log_and_print(f'###### DONE BUILDING CHARISData OBJECT FOR {self.object_name} #######')
        else:
            self.dataset = Dataset_PlaceHolder()
            append_dataset_info(build_charis_data, self.dataset)

        self.fake_fluxes = np.array(fake_fluxes)
        self.fake_seps = np.array(fake_seps)
        self.fake_fwhm = fake_fwhm
        self.fake_PAs = np.array(fake_PAs)

        self.trials = list()
        # START OF IF/ELSE AND FOR LOOPS FOR BUILDING TRIALS #
        if build_all_combos:
            if not isinstance(batched, tuple):
                for ani in annuli:
                    for subsec in subsections:
                        for mov in movement:
                            for spec in spectrum:
                                for cs in corr_smooth:
                                    for hp in highpass:
                                        self.trials.append(Trial(object_name=self.object_name, mask_xy=self.mask_xy,
                                                                 annuli=ani, subsections=subsec, movement=mov,
                                                                 numbasis=numbasis, spectrum=spec, corr_smooth=cs,
                                                                 fake_PAs=self.fake_PAs, fake_fluxes=self.fake_fluxes,
                                                                 fake_fwhm=self.fake_fwhm, fake_seps=self.fake_seps,
                                                                 rot_angs=self.dataset.PAs, flipx=self.dataset.flipx,
                                                                 dn_per_contrast=self.dataset.dn_per_contrast,
                                                                 wln_um=self.dataset.wvs, highpass=hp,
                                                                 length=self.dataset.leNgth))
            else:
                args = [annuli, subsections, movement, spectrum, corr_smooth, highpass]
                _, batchindex, batchsize = batched
                paramset = parameter_set_batcher(batchindex, batchsize, args)
                for params in paramset:
                    ani, subsec, mov, spec, cs, hp = params
                    self.trials.append(Trial(object_name=self.object_name, mask_xy=self.mask_xy,
                                             annuli=ani, subsections=subsec, movement=mov,
                                             numbasis=numbasis, spectrum=spec, corr_smooth=cs,
                                             fake_PAs=self.fake_PAs, fake_fluxes=self.fake_fluxes,
                                             fake_fwhm=self.fake_fwhm, fake_seps=self.fake_seps,
                                             rot_angs=self.dataset.PAs, flipx=self.dataset.flipx,
                                             dn_per_contrast=self.dataset.dn_per_contrast,
                                             wln_um=self.dataset.wvs, highpass=hp,
                                             length=self.dataset.leNgth))
        else:
            if not isinstance(batched, tuple):
                for ani, subsec, mov, spec, nb, cs, hp in zip(annuli, subsections, movement, spectrum,
                                                              numbasis, corr_smooth, highpass):
                    self.trials.append(Trial(object_name=self.object_name, mask_xy=self.mask_xy, annuli=ani,
                                             subsections=subsec, movement=mov, numbasis=[nb], spectrum=spec,
                                             corr_smooth=cs, fake_PAs=self.fake_PAs, fake_fluxes=self.fake_fluxes,
                                             fake_fwhm=self.fake_fwhm, fake_seps=self.fake_seps,
                                             rot_angs=self.dataset.PAs, flipx=self.dataset.flipx,
                                             dn_per_contrast=self.dataset.dn_per_contrast, wln_um=self.dataset.wvs,
                                             highpass=hp, length=self.dataset.leNgth))
            else:
                args = [annuli, subsections, movement, numbasis, spectrum, corr_smooth, highpass]
                _, batchindex, batchsize = batched
                startindex = batchsize * (batchindex - 1)  # 1-based indexing being passed in
                endindex = startindex + batchsize
                paramset = [arg[startindex: endindex] for arg in args]
                for params in zip(paramset[0], paramset[1], paramset[2], paramset[3], paramset[4], paramset[5],
                                  paramset[6]):
                    ani, subsec, mov, nb, spec, cs, hp = params
                    self.trials.append(Trial(object_name=self.object_name, mask_xy=self.mask_xy, annuli=ani,
                                             subsections=subsec, movement=mov, numbasis=[nb], spectrum=spec,
                                             corr_smooth=cs, fake_PAs=self.fake_PAs, fake_fluxes=self.fake_fluxes,
                                             fake_fwhm=self.fake_fwhm, fake_seps=self.fake_seps,
                                             rot_angs=self.dataset.PAs, flipx=self.dataset.flipx,
                                             dn_per_contrast=self.dataset.dn_per_contrast, wln_um=self.dataset.wvs,
                                             highpass=hp, length=self.dataset.leNgth))
        # END OF IF/ELSE AND FOR LOOPS #
        self.mode = mode
        self.overwrite = overwrite
        self.memorylite = memorylite
        self.write_to_log_and_print(f'############ DONE BUILDING TRIALS FOR {self.object_name} ############')

        if build_charis_data == 'temporary':  # removing from memory if just needed it for the rotation angles
            self.dataset = None

    def write_to_log(self, words, write_type='a'):
        with open(f'{self.object_name}/log.txt', write_type) as log_file:
            log_file.write(words)

    def write_to_log_and_print(self, words, write_type='a'):
        with open(f'{self.object_name}/log.txt', write_type) as log_file:
            log_file.write('\n' + words)
        print(words)

    def inject_fakes(self):
        if not isinstance(self.fake_fluxes[0], np.ndarray):  # regular like tutorial
            for fake_flux, sep in zip(self.fake_fluxes, self.fake_seps):
                flux_to_inject = fake_flux * self.dataset.dn_per_contrast  # UNcalibrating it
                for pa in self.fake_PAs:
                    inject_planet(frames=self.dataset.input, centers=self.dataset.centers,
                                  inputflux=flux_to_inject, astr_hdrs=self.dataset.wcs, radius=sep, pa=pa,
                                  fwhm=self.fake_fwhm)
        else:  # multiple tiers of planets
            for fluxgroup, pagroup in zip(self.fake_fluxes, self.fake_PAs):
                for fake_flux, sep in zip(fluxgroup, self.fake_seps):
                    flux_to_inject = fake_flux * self.dataset.dn_per_contrast  # UNcalibrating it
                    for pa in pagroup:
                        inject_planet(frames=self.dataset.input, centers=self.dataset.centers,
                                      inputflux=flux_to_inject, astr_hdrs=self.dataset.wcs, radius=sep,
                                      pa=pa, fwhm=self.fake_fwhm)
        self.write_to_log_and_print(f'############ DONE INJECTING FAKES FOR {self.object_name} ############')

    def run_KLIP_on_data_without_fakes(self, numthreads):
        if not os.path.exists(self.object_name):
            try:
                os.mkdir(self.object_name)
            except FileExistsError:
                # unnecessary 99% of the time, but once in a blue moon, I run into issues with this on Condor
                pass
        if not os.path.exists(self.object_name + '/klipped_cubes_Nfakes'):
            try:
                os.mkdir(self.object_name + '/klipped_cubes_Nfakes')
            except FileExistsError:
                # unnecessary 99% of the time, but once in a blue moon, I run into issues with this on Condor
                pass

        number_of_klip = len(self.trials)

        self.write_to_log_and_print('####### BEGINNING KLIP ON DATA WITHOUT FAKES #######\n'
                                    f'####### Number of KLIP Runs To Complete: {number_of_klip} #######\n')

        prevtime = time()
        for klip_runs, trial in enumerate(self.trials):  # klip_runs indicates how many have been previously completed
            # Update Every 20
            if klip_runs != 0 and klip_runs % 20 == 0:
                currenttime = time()
                minutes_per_run = round(((currenttime - prevtime) / 60) / 20, 2)
                self.write_to_log_and_print('####### {0}/{1} KLIP Runs Complete ({2}%) -- avg speed: {3} '
                                            'min/run #######'.format(klip_runs + 1, number_of_klip, round(float(
                    klip_runs + 1) / float(number_of_klip) * 100, 1), minutes_per_run))
                prevtime = time()

            if not self.overwrite:
                # not going to overwrite previous output
                filename = self.object_name + '/klipped_cubes_Nfakes' + self.object_name + '_withoutfakes_' + \
                           trial.klip_parameters + f'-KL{trial.numbasis}-speccube.fits'
                if os.path.exists(filename):
                    self.write_to_log_and_print(f"{filename} ALREADY EXISTS -- continuing without running KLIP on this "
                                                f"set of parameters")
                    continue

            with log_file_output(self.object_name):
                klip_dataset(self.dataset, outputdir=self.object_name + '/klipped_cubes_Nfakes',
                             fileprefix=self.object_name + '_withoutfakes_' + trial.klip_parameters,
                             annuli=trial.annuli, subsections=trial.subsections, movement=trial.movement,
                             numbasis=trial.numbasis, spectrum=trial.spectrum, corr_smooth=trial.corr_smooth,
                             highpass=trial.highpass, mode=self.mode, numthreads=numthreads, verbose=True,
                             lite=self.memorylite)

            # Update If Completely Done
            if (klip_runs + 1) == len(self.trials):
                self.write_to_log_and_print("\n### DONE WITH KLIP ON DATA WITH FAKES ###")

    def run_KLIP_on_data_with_fakes(self, numthreads):
        if not os.path.exists(self.object_name):
            try:
                os.mkdir(self.object_name)
            except FileExistsError:
                # unnecessary 99% of the time, but once in a blue moon, I run into issues with this on Condor
                pass
        if not os.path.exists(self.object_name + '/klipped_cubes_Wfakes'):
            try:
                os.mkdir(self.object_name + '/klipped_cubes_Wfakes')
            except FileExistsError:
                # unnecessary 99% of the time, but once in a blue moon, I run into issues with this on Condor
                pass

        number_of_klip = len(self.trials)

        self.write_to_log_and_print('####### BEGINNING KLIP ON DATA WITH FAKES #######\n'
                                    f'####### Number of KLIP Runs To Complete: {number_of_klip} #######\n')

        prevtime = time()
        for klip_runs, trial in enumerate(self.trials):  # klip_runs indicates how many have been previously completed
            if not self.overwrite:
                # not going to oveerwrite previous output
                filename = self.object_name + '/klipped_cubes_Wfakes' + self.object_name + '_withfakes_' + \
                           trial.klip_parameters + f'-KL{trial.numbasis}-speccube.fits'
                if os.path.exists(filename):
                    self.write_to_log_and_print(f"{filename} ALREADY EXISTS -- continuing without running KLIP on this "
                                                f"set of parameters")
                    continue

            with log_file_output(self.object_name):
                klip_dataset(self.dataset, outputdir=self.object_name + '/klipped_cubes_Wfakes',
                             fileprefix=self.object_name + '_withfakes_' + trial.klip_parameters,
                             annuli=trial.annuli, subsections=trial.subsections, movement=trial.movement,
                             numbasis=trial.numbasis, spectrum=trial.spectrum, corr_smooth=trial.corr_smooth,
                             highpass=trial.highpass, mode=self.mode, numthreads=numthreads, verbose=True,
                             lite=self.memorylite)

            # Update Every 20 or When Completely Done
            if klip_runs + 1 == len(self.trials):
                self.write_to_log_and_print("\n### DONE WITH KLIP ON DATA WITH FAKES ###")
            elif (klip_runs + 1) % 20 == 0:
                currenttime = time()
                minutes_per_run = round(((currenttime - prevtime) / 60) / 20, 2)
                self.write_to_log_and_print('####### {0}/{1} KLIP Runs Complete ({2}%) -- avg speed: {3} '
                                            'min/run #######'.format(klip_runs + 1, number_of_klip, round(float(
                    klip_runs + 1) / float(number_of_klip) * 100, 1), minutes_per_run))
                prevtime = time()

    def contrast_and_detection(self, run_contrast=True, run_planet_detection=True, datasetwithfakes=True):
        if not os.path.exists(self.object_name + '/detections') and run_planet_detection:
            try:
                os.mkdir(self.object_name + '/detections')
            except FileExistsError:
                # unnecessary 99% of the time, but once in a blue moon, I run into issues with this on Condor
                pass
        if not os.path.exists(self.object_name + '/calibrated_contrast') and datasetwithfakes:
            try:
                os.mkdir(self.object_name + '/calibrated_contrast')
            except FileExistsError:
                # unnecessary 99% of the time, but once in a blue moon, I run into issues with this on Condor
                pass
        if not os.path.exists(self.object_name + '/uncalibrated_contrast'):
            try:
                os.mkdir(self.object_name + '/uncalibrated_contrast')
            except FileExistsError:
                # unnecessary 99% of the time, but once in a blue moon, I run into issues with this on Condor
                pass

        if run_planet_detection and run_contrast:
            self.write_to_log_and_print(f"\n############## BEGINNING CONTRAST AND DETECTION FOR {self.object_name} "
                                        "##############")
        elif run_contrast:
            self.write_to_log_and_print(f"\n############## BEGINNING CONTRAST FOR {self.object_name} ##############")
        elif run_planet_detection:
            self.write_to_log_and_print(f"\n############## BEGINNING DETECTIONS FOR {self.object_name} ##############")

        if run_contrast:
            for i, trial in enumerate(self.trials):
                trial.get_contrast(contains_fakes=datasetwithfakes)
                if (i + 1) % 100 == 0:
                    print(f'# DONE WITH CONTRAST FOR {i + 1} TRIALS #')

        if run_planet_detection and run_contrast:
            self.write_to_log_and_print(f'### DONE WITH CONTRAST FOR {self.object_name}. BEGINNING DETECTION ###')
        elif run_contrast:
            self.write_to_log_and_print(f'### DONE WITH CONTRAST FOR {self.object_name}. ###')

        if run_planet_detection:
            # at the moment, can't parallelize because planet detection already utilizes (a little) parallelization
            for i, trial in enumerate(self.trials):
                trial.detect_planets(datasetwithfakes=datasetwithfakes, override=self.overwrite)
                if (i + 1) % 100 == 0:
                    print(f'# DONE WITH DETECTION FOR {i + 1} TRIALS #')

        if run_planet_detection:
            self.write_to_log_and_print(f"\n############## DONE WITH DETECTION FOR {self.object_name} ##############")