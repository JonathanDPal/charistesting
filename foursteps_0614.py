from astropy.io import fits
import numpy as np
from glob import glob
from pyklip.instruments.CHARIS import CHARISData
from copy import copy
from pyklip.parallelized import klip_dataset
from pyklip.klip import meas_contrast
import pylab as P
from csv import writer
from pyklip.fakes import retrieve_planet_flux
from pyklip.fakes import inject_planet
from astropy.wcs import WCS
from pyklip.kpp.utils.mathfunc import gauss2d
from pyklip.kpp.metric.crossCorr import calculate_cc
from pyklip.kpp.stat.statperPix_utils import get_image_stat_map_perPixMasking
from pyklip.kpp.detection.detection import point_source_detection
import pandas as pd
from math import ceil


def FWHMIOWA_calculator(speccubefile):
	"""
	Finds FWHM, IWA, and OWA for a opened CHARIS data cube.
	"""
	wavelengths = {'j': 1200e-9, 'h': 1550e-9, 'k': 2346e-9, 'broadband': 1550e-9}
	wavelength = wavelengths[str.lower(speccubefile[1].header['FILTNAME'])]
	D = 8
	lenslet_scale = 0.0162
	field_radius = 1.035
	FWHM = 2 * 1.22 * wavelength / D * 206265 / lenslet_scale
	IWA = 5
	OWA = (field_radius / lenslet_scale) - FWHM

	return FWHM, IWA, OWA


def calibrate_ss_contrast(speccubefile):
	"""
	Provides ratios to calibrate flux to be relative to the star at all wavelengths.
	"""
	cube = copy(speccubefile[1].data)

	# Number of Wavelength Bins
	Nwlm = cube.shape[0]

	# Creates Array of Wavelengths in nm
	lam_min = speccubefile[0].header['LAM_MIN']
	dloglam = speccubefile[0].header['DLOGLAM']
	wln_nm = lam_min * np.exp(np.arange(Nwlm)*dloglam)

	# Converts from micrometers to nanometers
	wln_um = wln_nm * 1.0e-3

	# Calculates the Spot/Star Ratio For Each Wavelength, in um
	spot_to_star = 2.72e-3 * (wln_um / 1.55) ** -2

	return wln_um, spot_to_star


class TestDataset:
	def __init__(self, fileset, object_name, mask_xy, fake_fluxes, fake_seps,
				 klipped_with_fakes_filepath, klipped_no_fakes_filepath, annuli,
				 subsections, movement, numbasis, mode='ADI+SDI', fake_fwhm=3.5, fake_PAs=[0,90,
																						   180,270]):
		# Creating CHARISData Object With UnKLIPped Data
		self.dataset_no_fakes = CHARISData(glob(fileset))
		self.dataset_with_fakes = copy(self.dataset_no_fakes)

		# Setting Object Name and Location
		self.object_name = object_name
		self.mask_xy = mask_xy

		# Info For Injecting (and later identifying) Fake Planets
		self.fake_fluxes = fake_fluxes
		self.fake_seps = fake_seps
		self.fake_fwhm = fake_fwhm
		self.fake_PAs = fake_PAs

		# KLIP Parameters
		self.annuli = annuli
		self.subsections = subsections
		self.movement = movement
		self.numbasis = numbasis
		self.mode = mode
		self.klip_parameters = [annuli, subsections, movement, numbasis, mode]

		# Where Contrast Info Will Be Stored
		self.uncalib_contrast = dict()
		self.calib_contrast = dict()
		self.contrast_seps = None

		# step two/three
		self.fakes_klipped_filepath = klipped_with_fakes_filepath
		self.no_fakes_klipped_filepath = klipped_no_fakes_filepath

		# step four
		self.detections_finders = None
		self.contrast_finders = None


	def injecting_fakes(self):
		"""
		Injects fake planets into CHARIS data.
		"""

		# Getting Values
		with fits.open(self.filelist[0]) as hdu:
			_, spot_to_star = calibrate_ss_contrast(hdu)

		# Inject Fake Planets
		for fake_flux, sep in zip(self.fake_fluxes, self.fake_seps):
			flux_to_inject = fake_flux * spot_to_star
			for pa in self.fake_PAs:
				inject_planet(self.dataset_with_fakes.input, self.dataset_with_fakes.centers,
							  flux_to_inject, self.dataset_with_fakes.wcs, sep, pa,
							  fwhm=self.fake_fwhm)


	def run_KLIP(self):
		# Running KLIP on Data Without Fakes
		klip_dataset(self.data_without_fakes, outputdir=self.no_fakes_klipped_filepath,
					 fileprefix=self.object_name+'_withoutfakes', annuli=self.annuli,
					 subsections=self.subsections, movement=self.movement, numbasis=self.numbasis,
					 mode=self.mode)

		# Running KLIP on Data With Fakes
		klip_dataset(self.data_with_fakes, outputdir=self.fakes_klipped_filepath,
					 fileprefix=object_name + '_withfakes', annuli=self.annuli,
					 subsections=self.subsections, movement=self.movement, numbasis=self.numbasis,
					 mode=self.mode)


	def ss_calibration(self, contains_fakes):
		"""
		Calibrates data cube with respect to star flux and provides a bunch of other information
		needed for contrast measurements.
		---
		Arg:
			contains_fakes (bool): Set to True if looking to calibrate data with the fakes and
									to False if looking to calibrate data without the fakes.

		Returns:
			wavelength in micrometers
			datacube calibrated with respect to star flux
			dataset center as a list of form [X-coor, Y-coor]
			dataset inner working angle
			dataset outer working angle
			dataset full width at half maximum
			dataset output WCS
		"""
		if contains_fakes:
			directory = self.fakes_klipped_filepath
			object_name_end = '_withfakes-KL'
		elif not contains_fakes:
			directory = self.no_fakes_klipped_filepath
			object_name_end = '_withoutfakes-KL'

		with fits.open(directory + '/' + self.object_name + object_name_end
					   + str(KL_modes_contrast_curve) + '-speccube.fits') as hdulist:
			wln_um, spot_to_star = calibrate_ss_contrast(hdulist)
			calib_cube = copy(hdulist[1].data) * spot_to_star[:, np.newaxis, np.newaxis]
			dataset_center = [hdulist[1].header['PSFCENTX'], hdulist[1].header['PSFCENTY']]
			dataset_fwhm, dataset_iwa, dataset_owa = FWHMIOWA_calculator(hdulist)
			output_wcs = WCS(header=hdulist[1].header, naxis=[1, 2])

		return wln_um, calib_cube, dataset_center, dataset_iwa, dataset_owa, dataset_fwhm, \
			   output_wcs


	def get_contrast(self, data_output_filepath, calib_cube, dataset_center, dataset_iwa,
					 dataset_owa, dataset_fwhm, output_wcs, wavelength_index=10,
					 contains_fakes=False):
		"""
		Measures contrast in an image and generates contrast curve.
		---
		Args:
			data_output_filepath (str): Filepath for contrast curve data.
										Should end in '.csv'
			calib_cube (np array): Data cube calibrated with regard to
									star/spot flux.
			dataset_center (list)
			dataset_iwa (float or int)
			dataset_owa (float)
			dataset_fwhm (float)
			output_wcs (instance of astropy.wcs WCS class)
			wavelength_index (int): Index of wavelength to use (subsets calibrated cube along
									wavelength axis).Default: 10
			contains_fakes (bool): Set to true if fakes present. Default: False
		"""

		frame = calib_cube[wavelength_index]

		if mask_xy is not None:
			x_pos = self.mask_xy[0]
			y_pos = self.mask_xy[1]

			ydat, xdat = np.indices(frame.shape)
			distance_from_planet = np.sqrt((xdat - x_pos) ** 2 + (ydat - y_pos) ** 2)
			frame[np.where(distance_from_planet <= 2 * dataset_fwhm)] = np.nan

		contrast_seps, contrast = meas_contrast(frame, dataset_iwa, dataset_owa, dataset_fwhm,
												center=dataset_center, low_pass_filter=True)

		# Calibrating For KLIP Subtraction If Fakes Present
		if contains_fakes:
			retrieved_fluxes = []
			for sep in injected_seps:
				fake_planet_fluxes = []
				for pa in injected_PAs:
					fake_flux = retrieve_planet_flux(frame, dataset_center, output_wcs, sep, pa,
													 searchrad=7)
					fake_planet_fluxes.append(fake_flux)
				retrieved_fluxes.append(np.mean(fake_planet_fluxes))

			algo_throughput = np.array(retrieved_fluxes) / np.array(injected_fluxes)

			correct_contrast = np.copy(contrast)
			for i, sep in enumerate(contrast_seps):
				closest_throughput_index = np.argmin(np.abs(sep - injected_seps))
				correct_contrast[i] /= algo_throughput[closest_throughput_index]

		# Saving Data
		if contains_fakes:
			self.calib_contrast = correct_contrast
			self.contrast_seps = contrast_seps
			with open(data_output_filepath, 'w+') as csvfile:
				csvwriter = writer(csvfile, delimiter=',')
				csvwriter.writerows([['Sep (Pixels)', 'Contrast']])
				csvwriter.writerows([contrast_seps, correct_contrast])
		elif not contains_fakes:
			with open(data_output_filepath, 'w+') as csvfile:
				csvwriter = writer(csvfile, delimiter=',')
				csvwriter.writerows([['Sep (Pixels)', 'Contrast']])
				csvwriter.writerows([contrast_seps, contrast])


	def step_three(self, SNR_threshold=3):
		"""
		Looks at a KLIPped dataset with fakes and indicates potential planets.
		---
		Args:
			SNR_threshold: Default: 3. In general, have this be the lowest value that we want to
							explore because it is super easy to just subset the output data to
							identify the subset that would have been identified at a higher SNR.
		"""

		with fits.open(self.fakes_klipped_filepath) as hdulist:
			image = hdulist[1].data
			center = [hdulist[0].header['PSFCENTX'], hdulist[0].header['PSFCENTY']]

		x_grid, y_grid = np.meshgrid(np.arange(-10,10), np.arange(-10,10))
		kernel_gauss = gauss2d(x_grid, y_grid)

		image_cc = calculate_cc(image, kernel_gauss, spectrum=np.ones(len(image[0])),
								nans2zero=True)

		SNR_map = get_image_stat_map_perPixMasking(image_cc, centroid=center, mask_radius=5,
												   Dr=2, type='SNR')

		detection_threshold = SNR_threshold
		pix2as = 1
		mask_radius = 15
		maskout_edge = 10

		candidates_table = point_source_detection(SNR_map, center, detection_threshold,
												  pix2as=pix2as, mask_radius=mask_radius,
												  maskout_edge=maskout_edge, IWA=None, OWA=None)

		# Write Out Detections To CSV File
		with open(output_prefix+'_SNR-'+str(SNR_threshold)+'.csv','w+') as csvfile:
			csvwriter = writer(csvfile, delimiter=',')
			csvwriter.writerows([['Index', 'SNR Value', 'PA', 'Sep (pix)',
								  'Sep (as)', 'x', 'y', 'row', 'col']])
			csvwriter.writerows(candidates_table)


	def step_four(self):
		"""
		Aggregates information and provides summary plots.
		"""

		detections_list = []
		for detection_finder in list(self.detections_finders):
			detections_list.append(glob(detection_finder))
		detections = dict()
		for i, filename in enumerate(detections_list):
			if filename[-6] == '-':
				try:
					snr_value = int(filename[-5])
				except ValueError:
					print("Filename {0} is not correctly formatted and will not be "
						  "included in the compiled data".format(filename))
					continue
			elif filename[-7] == '-':
				try:
					snr_value = int(filename[-6] + filename[-5])
				except ValueError:
					print("Filename {0} is not correctly formatted and will not be "
						  "included in the compiled data".format(filename))
					continue
			else:
				print("Filename {0} is not correctly formatted and will not be "
					  "included in the compiled data".format(filename))
				continue
			detections[i] = [pd.read_csv(filename), snr_value]

		snr_detections = dict()

		for i in list(detections.keys()):
			min_snr = detections[i][1]
			data = detections[i][0]
			max_snr = ceil(np.max(data['SNR Value']))
			snr_to_use = min_snr
			while snr_to_use <= max_snr:
				detections_at_snr_value = data[data['SNR Value'] > snr_to_use]
				for pa in detections_at_snr_value['PA']:


		contrasts_list = []
		for contrast_finder in list(self.contrast_finders):
			contrasts_list.append(glob(contrast_finder))
		contrasts = dict()
		for filename in contrasts_list:
			# Change the following section to fit some convention for saving
			# contrast data. Determine the convention and then include it in the
			# docstring for pt2 and probably even throw a test on it to make sure
			# that data_output_filename arg fits the convention.
			if filename[-6] == '-':
				whatever_value = filename[-5]
			elif filename[-7] == '-':
				whatever_value = filename[-6] + filename[-5]
			else:
				print("Filename {0} is not correctly formatted and will not be "
					  "included in the compiled data".format(filename))
				continue
			contrasts[whatever_value] = pd.read_csv(filename)



# Measuring Contrast, Generating Contrast Curve on Data Without Fakes &
# Saving the Contrast Data to a CSV File
no_fakes_wln_um, no_fakes_calib_cube, no_fakes_dataset_center, no_fakes_dataset_iwa, \
no_fakes_dataset_owa, no_fakes_dataset_fwhm, no_fakes_output_wcs = ss_calibration(False)

get_contrast(graph_output_filepath=filepath_uncalibrated_curve,
						data_output_filepath=filepath_uncal_contrast_data,
						wln_um=no_fakes_wln_um, calib_cube=no_fakes_calib_cube,
						dataset_center=no_fakes_dataset_center,
						dataset_iwa=no_fakes_dataset_iwa, dataset_owa=no_fakes_dataset_owa,
						dataset_fwhm=no_fakes_dataset_fwhm, output_wcs=no_fakes_output_wcs)

# Measuring Contrast, Generating Contrast Curve on Data With Fakes &
# Saving the Contrast Data to a CSV File
with_fakes_wln_um, with_fakes_calib_cube, with_fakes_dataset_center, \
with_fakes_dataset_iwa, with_fakes_dataset_owa, with_fakes_dataset_fwhm, \
with_fakes_output_wcs = ss_calibration(True)

get_contrast(graph_output_filepath=filepath_calibrated_curve,
						data_output_filepath=filepath_cal_contrast_data,
						wln_um=with_fakes_wln_um, calib_cube=with_fakes_calib_cube,
						dataset_center=with_fakes_dataset_center,
						dataset_iwa=with_fakes_dataset_iwa,
						dataset_owa=with_fakes_dataset_owa,
						dataset_fwhm=with_fakes_dataset_fwhm,
						output_wcs=with_fakes_output_wcs, contains_fakes=True,
						injected_fluxes=self.fake_fluxes, injected_seps=self.fake_seps,
						injected_PAs=self.fake_PAs)