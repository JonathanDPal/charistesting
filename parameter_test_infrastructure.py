from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from glob import glob
from pyklip.instruments.CHARIS import CHARISData
from copy import copy
from pyklip.parallelized import klip_dataset
from pyklip.klip import meas_contrast
import pylab as P
from csv import writer
from pyklip.fakes import inject_planet, retrieve_planet_flux
from pyklip.kpp.utils.mathfunc import gauss2d
from pyklip.kpp.metrics.crossCorr import calculate_cc
from pyklip.kpp.stat.statPerPix_utils import get_image_stat_map_perPixMasking
from pyklip.kpp.detection.detection import point_source_detection
import pandas as pd
from math import ceil
import os


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


# TestDataset Will Have a List of Trials Associated With It (one for each group of KLIP Parameters)
class Trial:
	def __init__(self, annuli, subsections, movement, numbasis, spectrum, mask_xy,
				 fake_PAs, fake_fluxes, object_name, fake_fwhm, fake_seps):
		self.annuli = annuli
		self.subsections = subsections
		self.movement = movement
		self.numbasis = numbasis
		self.spectrum = spectrum
		self.mask_xy = mask_xy
		self.fake_PAs = fake_PAs
		self.fake_fluxes = fake_fluxes
		self.object_name = object_name
		self.fake_fwhm = fake_fwhm
		self.fake_seps = fake_seps

		# String Identifying Parameters Used
		self.klip_parameters = str(annuli)+'Annuli_'+str(subsections)+'Subsections_'+str(
			movement)+'Movement_'+str(spectrum)+'Spectrum_'

		# Filepaths to KLIPped Datacubes
		self.filepaths_Wfakes = [self.object_name + '/klipped_cubes_Wfakes/' + self.object_name + \
									'_withfakes_' + self.klip_parameters +'-KL{0}-' \
									'speccube.fits'.format(nb) for nb in self.numbasis]
		self.filepaths_Nfakes = [self.object_name + '/klipped_cubes_Nfakes/' + self.object_name + \
									'_withoutfakes_' + self.klip_parameters +'-KL{0}-' \
									'speccube.fits'.format(nb) for nb in self.numbasis]

		# Filepath to Save
		self.filepath_detections_prefixes = [self.object_name + '/detections/{0}_KL{1}_SNR-'.format(
			self.klip_parameters, nb) for nb in self.numbasis]
		# Setting Up Filepath
		if not os.path.exists(self.object_name):
			os.mkdir(self.object_name)
		if not os.path.exists(self.object_name + '/detections'):
			os.mkdir(self.object_name + '/detections')

		# Data To Be Added Later
		self.calib_cube = None
		self.uncalib_contrast = None
		self.calib_contrast = None
		self.detections = None
		self.classified_detections = None


	def get_contrast(self, contains_fakes, wavelength_index=10):
		"""
		Measures contrast in an image; saves contrast data both to a CSV file and to the object.
		---
		Args:
			wavelength_index (int): Index of wavelength to use (subsets calibrated cube along
									wavelength axis).Default: 10
			contains_fakes (bool): Set to true if fakes present. Default: False
		"""
		if contains_fakes:
			filepaths = self.filepaths_Wfakes
		else:
			filepaths = self.filepaths_Nfakes

		for i, filepath in enumerate(filepaths):
			with fits.open(filepath) as hdulist:
				wln_um, spot_to_star = calibrate_ss_contrast(hdulist)
				calib_cube = copy(hdulist[1].data) * spot_to_star[:, np.newaxis, np.newaxis]
				dataset_center = [hdulist[1].header['PSFCENTX'], hdulist[1].header['PSFCENTY']]
				dataset_fwhm, dataset_iwa, dataset_owa = FWHMIOWA_calculator(hdulist)
				output_wcs = WCS(header=hdulist[0].header, naxis=[1, 2])

			frame = calib_cube[wavelength_index]

			if self.mask_xy is not None:
				x_pos = self.mask_xy[0]
				y_pos = self.mask_xy[1]

				ydat, xdat = np.indices(frame.shape)
				distance_from_planet = np.sqrt((xdat - x_pos) ** 2 + (ydat - y_pos) ** 2)
				frame[np.where(distance_from_planet <= 2 * dataset_fwhm)] = np.nan

			contrast_seps, contrast = meas_contrast(frame, dataset_iwa, dataset_owa,
													dataset_fwhm, center=dataset_center,
													low_pass_filter=True)

			# Calibrating For KLIP Subtraction If Fakes Present
			if contains_fakes:
				retrieved_fluxes = []
				for sep in self.fake_seps:
					fake_planet_fluxes = []
					for pa in self.fake_PAs:
						fake_flux = retrieve_planet_flux(frame, dataset_center,
														 output_wcs, sep, pa,
														 searchrad=7)
						fake_planet_fluxes.append(fake_flux)
					retrieved_fluxes.append(np.mean(fake_planet_fluxes))

				algo_throughput = np.array(retrieved_fluxes) / np.array(self.fake_fluxes)

				correct_contrast = np.copy(contrast)
				for j, sep in enumerate(contrast_seps):
					closest_throughput_index = np.argmin(np.abs(sep - self.fake_seps))
					correct_contrast[j] /= algo_throughput[closest_throughput_index]

			# Making Sure That Directories Exist For Saving Data
			if not os.path.exists(self.object_name):
				os.mkdir(self.object_name)
			if not os.path.exists(self.object_name + '/calibrated_contrast'):
				os.mkdir(self.object_name + '/calibrated_contrast')
			if not os.path.exists(self.object_name + '/uncalibrated_contrast'):
				os.mkdir(self.object_name + '/uncalibrated_contrast')

			# Saving Data to Object and to CSV File
			if contains_fakes:
				self.calib_contrast = [contrast_seps, correct_contrast]
				data_output_filepath = self.object_name + \
									   '/calibrated_contrast/{0}_KL{1}_contrast.csv'.format(
										   self.klip_parameters, self.numbasis[i])
				with open(data_output_filepath, 'w+') as csvfile:
					csvwriter = writer(csvfile, delimiter=',')
					csvwriter.writerows([['Sep (Pixels)', 'Contrast']])
					csvwriter.writerows([contrast_seps, correct_contrast])
			else:
				self.uncalib_contrast = [contrast_seps, contrast]
				data_output_filepath = self.object_name + \
									   '/uncalibrated_contrast/{0}_KL{1}_contrast.csv'.format(
										   self.klip_parameters, self.numbasis[i])
				with open(data_output_filepath, 'w+') as csvfile:
					csvwriter = writer(csvfile, delimiter=',')
					csvwriter.writerows([['Sep (Pixels)', 'Contrast']])
					csvwriter.writerows([contrast_seps, contrast])


	def detect_planets(self, SNR_threshold=3):
		"""
		Looks at a KLIPped dataset with fakes and indicates potential planets.
		---
		Args:
			SNR_threshold: Default: 3. In general, have this be the lowest value that we want to
							explore because it is super easy to just subset the output data to
							identify the subset that would have been identified at a higher SNR.
		"""
		for i, filepath in enumerate(self.filepaths_Wfakes):
			with fits.open(filepath) as hdulist:
				image = hdulist[1].data
				center = [hdulist[1].header['PSFCENTX'], hdulist[1].header['PSFCENTY']]

			x_grid, y_grid = np.meshgrid(np.arange(-10,10), np.arange(-10,10))
			kernel_gauss = gauss2d(x_grid, y_grid)

			# flat spectrum given so that it collapses it into one image, instead of giving seperate
			# images for each wavelength
			image_cc = calculate_cc(image, kernel_gauss, spectrum=np.ones(len(image[0])),
									nans2zero=True)

			SNR_map = get_image_stat_map_perPixMasking(image_cc, centroid=center, mask_radius=5,
													   Dr=2, type='SNR')

			candidates_table = point_source_detection(SNR_map, center, SNR_threshold,
													  pix2as=1, mask_radius=15,
													  maskout_edge=10, IWA=None, OWA=None)

			self.detections = candidates_table

			candidates = pd.DataFrame(candidates_table, columns=['Index', 'SNR Value', 'PA',
																'Sep (pix)', 'Sep (as)', 'x',
																'y', 'row', 'col'])
			real_planet = []
			for _, row in candidates.iterrows():
				if np.min(np.array(row['PA']) - np.array(
						self.fake_PAs)) > 0.5 * self.fake_fwhm or np.min(
						np.array(row['Sep (pix)']) - np.array(self.fake_seps)) > 2 or np.min(
						np.array(row['Sep (as)']) - np.array(self.fake_seps)) > 2:
						real_planet.append(False)
				else:
						real_planet.append(True)
			candidates['Injected'] = real_planet
			self.classified_detections = candidates
			candidates.to_csv('{0}{1}.csv'.format(self.filepath_detections_prefixes[i],
												   str(SNR_threshold)))


	def __eq__(self, other):
		"""
		Checks to see if two Trials have the same KLIP parameters. Intended for testing
		out code functionality.
		"""
		return self.annuli == other.annuli and self.subsections == other.subsections and \
			   self.movement == other.movement and self.numbasis == other.numbasis and \
			   self.spectrum == other.spectrum and self.mode == other.mode


# Each Object (eg. HD1160, BetaPic) Will Have An Instance of TestDataset Associated With It
class TestDataset:
	def __init__(self, fileset, object_name, mask_xy, fake_fluxes, fake_seps,
				 annuli, subsections, movement, numbasis, spectrum=['methane', None],
				 fake_fwhm=3.5, fake_PAs=[0,90,180,270]):
		"""
		Args:
			fileset: Something probably going like 'directory/*.fits' to let glob find files.
			object_name: String
			mask_xy: [X-coor, Y-coor]
			fake_fluxes: Integer or List of Integers. Must be same length as fake_seps.
			fake_seps: Integer or List of Integers. Must be same length as fake_fluxes.
			annuli: Integer or List of Integers
			subsections: Integer or List of Integers
			movement: Integer or List of Integers
			numbasis: Integer or List of Integers
			spectrum: Either 'methane' or None
			fake_fwhm: The FWHM for the injected PSF for fake planets
			fake_PAs: Integer or List of Integers.
		"""
		# Setting Object Name and Location
		self.object_name = object_name
		self.mask_xy = mask_xy

		print("##################### STARTING WORK ON {0} #####################".format(self.object_name))

		# Creating CHARISData Object With UnKLIPped Data
		self.fileset = glob(fileset)
		self.dataset_no_fakes = CHARISData(self.fileset)
		print("####### DONE BUILDING CHARISData OBJECT FOR {0} ########".format(self.object_name))
		self.dataset_with_fakes = copy(self.dataset_no_fakes)

		# Info For Injecting (and later identifying) Fake Planets
		self.fake_fluxes = fake_fluxes
		self.fake_seps = fake_seps
		self.fake_fwhm = fake_fwhm
		self.fake_PAs = fake_PAs

		# Building Trials
		self.trials = []
		for ani in list(annuli):
			for subsec in list(subsections):
				for mov in list(movement):
						for spec in list(spectrum):
								self.trials.append(Trial(annuli=ani, subsections=subsec,
														 movement=mov, numbasis=numbasis,
														 spectrum=spec,
														 mask_xy=mask_xy, fake_PAs=fake_PAs,
														 fake_fluxes=fake_fluxes,
														 object_name=object_name,
														 fake_fwhm=fake_fwhm, fake_seps=fake_seps))
		print("############## DONE BUILDING TRIALS FOR {0} ##############".format(self.object_name))


	def inject_fakes(self):
		"""
		Injects fake planets into CHARIS data.
		"""

		# Getting Values
		with fits.open(self.fileset[0]) as hdu:
			_, spot_to_star = calibrate_ss_contrast(hdu)

		# Inject Fake Planets
		for fake_flux, sep in zip(self.fake_fluxes, self.fake_seps):
			flux_to_inject = fake_flux * spot_to_star
			for pa in self.fake_PAs:
				inject_planet(self.dataset_with_fakes.input, self.dataset_with_fakes.centers,
							  flux_to_inject, self.dataset_with_fakes.wcs, sep, pa,
							  fwhm=self.fake_fwhm)

		print("############## DONE INJECTING FAKES FOR {0} ##############".format(self.object_name))


	def run_KLIP(self):
		# Making Sure Output Directories Exist
		if not os.path.exists(self.object_name):
			os.mkdir(self.object_name)
		if not os.path.exists(self.object_name+'/klipped_cubes_Wfakes'):
			os.mkdir(self.object_name+'/klipped_cubes_Wfakes')
		if not os.path.exists(self.object_name+'/klipped_cubes_Nfakes'):
			os.mkdir(self.object_name+'/klipped_cubes_Nfakes')

		print("############## BEGINNING KLIP FOR {0} ##############".format(self.object_name))
		print("####### Total KLIP Runs to Complete: {0} #######".format(len(self.trials) * 2))

		for i, trial in enumerate(self.trials):
			# Running KLIP on Data With Fakes
			klip_dataset(self.dataset_with_fakes, outputdir=self.object_name+'/klipped_cubes_Wfakes',
						 fileprefix=self.object_name + '_withfakes_' + trial.klip_parameters,
						 annuli=trial.annuli,subsections=trial.subsections,
						 movement=trial.movement, numbasis=trial.numbasis, spectrum=trial.spectrum, verbose=False)

			# Running KLIP on Data Without Fakes
			klip_dataset(self.dataset_no_fakes, outputdir=self.object_name+'/klipped_cubes_Nfakes',
						 fileprefix=self.object_name+'_withoutfakes_' + trial.klip_parameters,
						 annuli=trial.annuli, subsections=trial.subsections,
						 movement=trial.movement, numbasis=trial.numbasis, spectrum=trial.spectrum, verbose=False)
			if i + 1 == len(self.trials):
				print("############## DONE WITH KLIP FOR {0} ##############".format(self.object_name))
			elif ((i+1) * 2) % 10 == 0:
				print("####### {0}/{1} KLIP Runs Complete ({2}%) #######".format((i+1) * 2,
						len(self.trials) * 2, round(float(i+1)/float(len(self.trials)), 3) * 100))



	def contrast_and_detection(self):
		print("############## BEGINNING CONTRAST AND DETECTION FOR {0} ##############".format(
			self.object_name))
		for i, trial in enumerate(self.trials):
		# 	for tf in [True, False]:
		# 		trial.get_contrast(tf)
			trial.detect_planets()
			if i + 1 == len(self.trials):
				print('############## DONE WITH CONTRAST AND DETECTION FOR {0} '
					  '##############'.format(self.object_name))
			elif ((i+1) * 2) % 10 == 0:
				print("####### Detection and contrast complete for {0}/{1} Trials ({2}%) "
					  "#######".format((i+1) * 2, len(self.trials) * 2, round(float(i+1)/float(len(
					self.trials)), 3) * 100))
