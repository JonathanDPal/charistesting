from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from glob import glob
from pyklip.instruments.CHARIS import CHARISData
from copy import copy, deepcopy
from pyklip.parallelized import klip_dataset
from pyklip.klip import meas_contrast
from csv import writer
from pyklip.fakes import inject_planet, retrieve_planet_flux
from pyklip.kpp.utils.mathfunc import gauss2d
from pyklip.kpp.metrics.crossCorr import calculate_cc
from pyklip.kpp.stat.statPerPix_utils import get_image_stat_map_perPixMasking
from pyklip.kpp.detection.detection import point_source_detection
import pandas as pd
import sys, os
import matplotlib.pyplot as plt
from contextlib import contextmanager
import inspect
import warnings


@contextmanager
def log_file_output(directory, write_type='a'):
	"""
	Has outputs written out to a log file in the specified directory instead of printed in terminal.
	"""
	with open('{0}/log.txt'.format(directory), write_type) as log_file:
		old_stdout = sys.stdout
		sys.stdout = log_file
		try:
			yield
		finally:
			sys.stdout = old_stdout


def FWHMIOWA_calculator(speccubefile, filtname=None):
	"""
	Finds FWHM, IWA, and OWA for a opened CHARIS data cube. Thanks to Dr. Tobin for this.
	(https://docs.google.com/document/d/1S1Oo9QweKwnOfmv6fu28bb75lYeQXGzn/edit)
	"""
	wavelengths = {'j': 1200e-9, 'h': 1550e-9, 'k': 2346e-9, 'broadband': 1550e-9}
	if filtname is None:
		wavelength = wavelengths[str.lower(speccubefile[1].header['FILTNAME'])]
	else:
		wavelength = wavelengths[str.lower(filtname)]
	D = 8
	lenslet_scale = 0.0162
	field_radius = 1.035
	FWHM = 2 * 1.22 * wavelength / D * 206265 / lenslet_scale
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

	# Returns modified dataset object
	return dataset


def pasep_to_xy(PAs, seps):
	"""
	Takes lists of position angles and seperations and yields a numpy array with x-y coordinates for each combo, using
	the convention used in the table of values outputted by the planet detection software. The origin used in their
	convention is the center of the star's PSF.
	"""
	PAs = [float(pa) for pa in PAs]
	seps = [float(sep) for sep in seps]
	radians = np.array(PAs) / 180 * np.pi
	locs = []
	for sep in seps:
		for rad in radians:
			x = -np.sin(rad) * sep
			y = np.cos(rad) * sep
			loc = [x, y]
			locs.append(loc)
	return locs


# TestDataset Will Have a List of Trials Associated With It (one for each group of KLIP Parameters)
class Trial:
	"""
	NOTE: The user will almost certainly not interact with this class directly, rather they will interact with an
	instance of TestDataset and that instance of TestDataset will interact with instances of this class.
	---
	Stores a particular set of KLIP parameters and then is able to run contrast measurement or planet detection code
	from KLIP for the KLIP output with that particular set of parameters.
	"""
	def __init__(self, object_name, mask_xy, annuli, subsections, movement, numbasis, spectrum, corr_smooth,
				 fake_PAs, fake_fluxes, fake_fwhm, fake_seps, dn_per_contrast, wln_um, highpass, length):
		self.object_name = object_name
		self.mask_xy = mask_xy

		self.annuli = annuli
		self.subsections = subsections
		self.movement = movement
		self.numbasis = numbasis
		self.spectrum = spectrum
		self.corr_smooth = corr_smooth

		self.fake_PAs = np.array(fake_PAs)
		self.fake_fluxes = fake_fluxes
		self.fake_fwhm = fake_fwhm
		self.fake_seps = np.array(fake_seps)

		self.dn_per_contrast = dn_per_contrast
		self.wln_um = wln_um

		# Switching Highpass To Original Form If Necessary
		if isinstance(highpass, (int, float)) and not isinstance(highpass, bool):
			highpass = float(highpass)
			self.highpass = length / (highpass * 2 * np.sqrt(2 * np.log(2)))
		else:
			self.highpass = highpass

		# This Value Will Be Written In After KLIP Run Performed
		self.output_wcs = None

		# String Identifying Parameters Used (Used Later For Saving Contrast Info)
		self.klip_parameters = str(annuli)+'Annuli_'+str(subsections)+'Subsections_'+str(movement)+'Movement_'+str(
								spectrum)+'Spectrum_'+str(corr_smooth)+'Smooth_'+str(highpass)+'Highpass_'

		# Filepaths to KLIPped Datacubes
		self.filepaths_Wfakes = [self.object_name + '/klipped_cubes_Wfakes/' + self.object_name + '_withfakes_' +
								 self.klip_parameters +'-KL{0}-speccube.fits'.format(nb) for nb in self.numbasis]
		self.filepaths_Nfakes = [self.object_name + '/klipped_cubes_Nfakes/' + self.object_name + '_withoutfakes_' +
								 self.klip_parameters +'-KL{0}-speccube.fits'.format(nb) for nb in self.numbasis]

		# Filepath to Save Planet Detection Output To
		self.filepath_detections_prefixes = [self.object_name + '/detections/{0}_KL{1}_SNR-'.format(
												self.klip_parameters, nb) for nb in self.numbasis]


	def get_contrast(self, contains_fakes, wavelength_index=10):
		"""
		Measures contrast, then saves a contrast curve as a PNG and the contrast data as a CSV.
		---
		Args:
			contains_fakes (bool): Set to True if looking to make measurements on data with fakes, vice versa.
			wavelength_index (int): Index of wavelength to look at for contrast information. Default: 1.63 um
		"""
		if contains_fakes:
			filepaths = self.filepaths_Wfakes
		else:
			filepaths = self.filepaths_Nfakes

		# Measuring Contrast For Each Set of KL Modes
		for filepath_index, filepath in enumerate(filepaths): # filepath_index used to identify number of KL modes
			with fits.open(filepath) as hdulist:
				cube = copy(hdulist[1].data)
				dataset_center = [hdulist[1].header['PSFCENTX'], hdulist[1].header['PSFCENTY']]
				dataset_fwhm, dataset_iwa, dataset_owa = FWHMIOWA_calculator(hdulist)

			# Taking Slice of Cube and then Calibrating It
			frame = cube[wavelength_index]
			frame /= self.dn_per_contrast[wavelength_index]

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

			# Applying Mask to Fake Planets
			if contains_fakes:
				fakelocs = pasep_to_xy(self.fake_PAs, self.fake_seps)
				for fl in fakelocs:
					x_pos = fl[0]
					y_pos = fl[1]

					ydat, xdat = np.indices(frame.shape)
					distance_from_planet = np.sqrt((xdat - x_pos) ** 2 + (ydat - y_pos) ** 2)
					frame[np.where(distance_from_planet <= 2 * dataset_fwhm)] = np.nan

			# Measuring Contrast
			contrast_seps, contrast = meas_contrast(frame, dataset_iwa, dataset_owa, dataset_fwhm,
													center=dataset_center, low_pass_filter=True)

			# Calibrating For KLIP Subtraction If Fakes Present
			if contains_fakes:
				retrieved_fluxes = []
				for sep in self.fake_seps:
					fake_planet_fluxes = []
					for pa in self.fake_PAs:
						fake_flux = retrieve_planet_flux(frame, dataset_center, self.output_wcs, sep, pa, searchrad=7)
						fake_planet_fluxes.append(fake_flux)
					retrieved_fluxes.append(np.mean(fake_planet_fluxes))

				numgroups = len(self.fake_seps)
				groupsize = int(len(self.fake_fluxes) / len(self.fake_seps))
				fluxgroups = [[self.fake_fluxes[i * (groupsize + 1) + j] for i in range(groupsize)] for j in range(
					numgroups)]
				fluxes = []
				for i in range(numgroups):
					fxs = []
					for j in range(groupsize):
						fxs.append(fluxgroups[i][j])
					fluxes.append(np.mean(fxs))
				algo_throughput = np.array(retrieved_fluxes) / np.array(fluxes)

				correct_contrast = np.copy(contrast)
				for j, sep in enumerate(contrast_seps):
					closest_throughput_index = np.argmin(np.abs(sep - self.fake_seps))
					correct_contrast[j] /= algo_throughput[closest_throughput_index]

			# Always Going to Save Uncalibrated Contrast
			if not os.path.exists(self.object_name + '/uncalibrated_contrast'):
				os.mkdir(self.object_name + '/uncalibrated_contrast')
			data_output_filepath = self.object_name + '/uncalibrated_contrast/{0}_KL{1}_contrast.csv'.format(
				self.klip_parameters, self.numbasis[filepath_index])
			df = pd.DataFrame()
			df['Seperation'] = contrast_seps
			df['Uncalibrated Contrast'] = contrast
			wavelength = round(self.wln_um[wavelength_index], 2)
			title = 'Uncalibrated Contrast at {0}um ({1})'.format(wavelength, self.object_name)
			df.plot(x='Seperation', y='Uncalibrated Contrast', legend=False, title=title)
			plt.ylabel('Uncalibrated Contrast')
			plt.xlabel('Seperation')
			plt.semilogy()
			plt.savefig(data_output_filepath[0:-4] + '.png')
			df.to_csv(data_output_filepath)

			# If Fakes Present, Then Use Them To Calibrate Contrast
			if contains_fakes:
				if not os.path.exists(self.object_name + '/calibrated_contrast'):
					os.mkdir(self.object_name + '/calibrated_contrast')
				data_output_filepath = self.object_name + '/calibrated_contrast/{0}_KL{1}_contrast.csv'.format(
										   self.klip_parameters, self.numbasis[filepath_index])
				df = pd.DataFrame()
				df['Seperation'] = contrast_seps
				df['Calibrated Contrast'] = contrast
				wavelength = round(self.wln_um[wavelength_index], 2)
				title = 'Calibrated Contrast at {0}um'.format(wavelength)
				df.plot(x='Seperation', y='Calibrated Contrast', legend=False, title=title)
				plt.ylabel('Calibrated Contrast')
				plt.xlabel('Seperation')
				plt.semilogy()
				plt.savefig(data_output_filepath[0:-4] + '.png')
				df.to_csv(data_output_filepath)


	def detect_planets(self, SNR_threshold=2, datasetwithfakes=True):
		"""
		Looks at a KLIPped dataset with fakes and indicates potential planets.
		---
		Args:
			SNR_threshold: Default: 2. Set this to the lowest value to be looked at. Later on, the subsets of
										detections can be created for other SNR_thresholds.
			datasetwithfakes (Bool): Default: True. If True, then run planet detection on dataset containing injected
										planets; if False, then run planet detection on dataset not containing
										injected planets.
		"""
		if datasetwithfakes:
			filepaths = self.filepaths_Wfakes
		else:
			filepaths = self.filepaths_Nfakes

		for filepath_index, filepath in enumerate(filepaths): # filepath_index used to identify number of KL modes
			with fits.open(filepath) as hdulist:
				image = copy(hdulist[1].data)
				center = [hdulist[1].header['PSFCENTX'], hdulist[1].header['PSFCENTY']]

			x_grid, y_grid = np.meshgrid(np.arange(-10,10), np.arange(-10,10))
			kernel_gauss = gauss2d(x_grid, y_grid)

			# flat spectrum given here for generating cross-coorelated image so that pyKLIP collapses it into one
			# image, instead of giving seperate images for each wavelength
			image_cc = calculate_cc(image, kernel_gauss, spectrum=np.ones(len(image[0])), nans2zero=True)

			SNR_map = get_image_stat_map_perPixMasking(image_cc, centroid=center, mask_radius=5, Dr=2, type='SNR')

			candidates_table = point_source_detection(SNR_map, center, SNR_threshold, pix2as=1, mask_radius=15,
													  maskout_edge=10, IWA=None, OWA=None)

			candidates = pd.DataFrame(candidates_table, columns=['Index', 'SNR Value', 'PA', 'Sep (pix)',
																 'Sep (as)', 'x', 'y', 'row', 'col'])
			injected = [] # going to be an additional column of candidates DataFrame
			fakelocs = pasep_to_xy(self.fake_PAs, self.fake_seps) # where planets were injected

			# Determining Whether Candidates Are Planets or Not
			xpos = candidates['x']
			ypos = candidates['y']''
			candidate_locations = zip(xpos, ypos)

			def distance(xy1, xy2):
				return np.sqrt((xy1[0] - xy2[0]) ** 2 + (xy1[1] - xy2[1]) ** 2)

			if not isinstance(self.mask_xy[0], (list, tuple)):
				self.mask_xy = [self.mask_xy]

			distances_from_fakes = []
			distances_from_targets = []
			for c in candidate_locations:
				distances = []
				for fl in fakelocs:
					distances.append(distance(c, fl))
				distances_from_fakes.append(np.min(distances))
				distances2 = []
				for mask in self.mask_xy:
					mask = np.array(mask) - 100 # aligning coordinate systems
					distances2.append(distance(c, mask))
				distances_from_targets.append(np.min(distances2))

			for d1, d2 in zip(distances_from_targets, distances_from_fakes):
				# If Detection Within 1 FWHM of Location, Considered Legit
				if d1 < self.fake_fwhm:
					injected.append("Science Target")
				elif d2 < self.fake_fwhm:
					injected.append(True)
				else:
					injected.append(False)

			# Appending Information to Previous candidates DataTable
			candidates['Distance From Fakes'] = distances_from_fakes
			candidates['Distance From Targets'] = distances_from_targets
			candidates['Injected'] = injected

			# Saving Information
			candidates.to_csv('{0}{1}.csv'.format(self.filepath_detections_prefixes[filepath_index],
												  str(SNR_threshold)))


	def __eq__(self, other):
		"""
		Checks to see if two Trials have the same KLIP parameters. Intended for testing out code functionality.
		Thanks to akritigoswami for some of this code.
		(https://www.geeksforgeeks.org/how-to-get-a-list-of-class-attributes-in-python/)
		"""
		equal_attributes = list()
		for i, j in zip(inspect.getmembers(self), inspect.getmembers(other)):
			if i[0].startswith('_') or inspect.ismethod(i[1]):
				continue
			else:
				equal_attributes.append(i[1] == j[1])
				if i[1] != j[1]:
					break
		return np.sum(equal_attributes) == len(equal_attributes)


# Each Object (eg. HD1160, BetaPic) Will Have An Instance of TestDataset Associated With It
class TestDataset:
	"""
	The main object which the user will interact with. Will load in CHARIS fileset into CHARISData class (see
	pyklip.instruments.CHARIS) and then create an instance of Trial for each set of KLIP parameters to be looked at.
	"""
	def __init__(self, fileset, object_name, mask_xy, fake_fluxes, fake_seps, annuli, subsections, movement,
				 numbasis, corr_smooth, highpass, spectrum, fake_fwhm, fake_PAs, mode):
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
			corr_smooth:
			highpass:
			spectrum: Either 'methane' or None_
			fake_fwhm: The FWHM for the injected PSF for fake planets
			fake_PAs: Integer or List of Integers.
		"""
		# Setting Object Name and Location
		self.object_name = object_name
		self.mask_xy = mask_xy

		# Making Necessary Directory
		if not os.path.exists(self.object_name):
			os.mkdir(self.object_name)

		# Getting Started On Log File
		self.write_to_log('Title for Set: {0}'.format(object_name), 'w')
		self.write_to_log('\nFileset: {0}'.format(fileset))
		param_names = ['Annuli', 'Subsections', 'Movement', 'Numbasis', 'Corr_Smooth', 'Highpass', 'Spectrum',
					   'Mode', 'Fake Fluxes', 'Fake Seps', 'Fake PAs', 'Fake FWHM']
		params = [annuli, subsections, movement, numbasis, spectrum, corr_smooth, highpass]
		number_of_trials = np.prod([len(p) for p in params])
		params.append(mode)
		params.append(fake_fluxes)
		params.append(fake_seps)
		params.append(fake_PAs)
		params.append(fake_fwhm)
		self.write_to_log('\nNumber of Trials: {0}'.format(number_of_trials))
		for name, param in zip(param_names, params):
			self.write_to_log('\n{0}: {1}'.format(name, param))

		self.write_to_log_and_print("############### STARTING WORK ON {0} ################\n".format(self.object_name))

		# Creating CHARISData Object With UnKLIPped Data
		self.fileset = glob(fileset)
		with log_file_output(self.object_name):
			self.dataset = make_dn_per_contrast(CHARISData(self.fileset))

		self.write_to_log_and_print("###### DONE BUILDING CHARISData OBJECT FOR {0} #######".format(self.object_name))

		self.length = self.dataset.input.shape[1]

		# Info For Injecting (and later identifying) Fake Planets
		self.fake_fluxes = fake_fluxes
		self.fake_seps = fake_seps
		self.fake_fwhm = fake_fwhm
		self.fake_PAs = fake_PAs

		# Building Trials
		self.trials = []
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
															 dn_per_contrast=self.dataset.dn_per_contrast,
															 wln_um=self.dataset.wvs, highpass=hp,
															 length=self.length))

		self.mode = mode
		self.write_to_log_and_print("############ DONE BUILDING TRIALS FOR {0} ############".format(self.object_name))


	def write_to_log(self, words, write_type='a'):
		with open('{0}/log.txt'.format(self.object_name), write_type) as log_file:
			log_file.write(words)


	def write_to_log_and_print(self, words, write_type='a'):
		with open('{0}/log.txt'.format(self.object_name), write_type) as log_file:
			log_file.write('\n'+words)
		print(words)


	def inject_fakes(self):
		# Inject Fake Planets
		if len(self.fake_fluxes) == len(self.fake_seps):
			for fake_flux, sep in zip(self.fake_fluxes, self.fake_seps):
				flux_to_inject = fake_flux * self.dataset.dn_per_contrast # UNcalibrating it
				for pa in self.fake_PAs:
					inject_planet(frames=self.dataset.input, centers=self.dataset.centers,
								  inputflux=flux_to_inject, astr_hdrs=self.dataset.wcs, radius=sep, pa=pa,
								  fwhm=self.fake_fwhm)
		elif len(self.fake_fluxes) % len(self.fake_seps) == 0:
			numgroups = len(self.fake_seps)
			groupsize = int(len(self.fake_fluxes) / len(self.fake_seps))
			fluxes = [[self.fake_fluxes[i * (groupsize + 1) + j] for i in range(groupsize)] for j in range(
				numgroups)]
			pas = [[self.fake_PAs[i * (groupsize + 1) + j] for i in range(groupsize)] for j in range(
				numgroups)]
			for i in range(numgroups):
				for fake_flux, sep in zip(fluxes[i], self.fake_seps):
					flux_to_inject = fake_flux * self.dataset.dn_per_contrast  # UNcalibrating it
					for pa in pas[i]:
						inject_planet(frames=self.dataset.input, centers=self.dataset.centers,
									  inputflux=flux_to_inject, astr_hdrs=self.dataset.wcs, radius=sep,
									  pa=pa, fwhm=self.fake_fwhm)
		else:
			self.write_to_log_and_print("inject_fakes function called, but fake_fluxes and fake_seps are not "
										"compatible, so no injection occurred.")

		self.write_to_log_and_print("############ DONE INJECTING FAKES FOR {0} ############".format(self.object_name))


	def run_KLIP_on_data_without_fakes(self):
		# Making Sure Output Directories Exist
		if not os.path.exists(self.object_name):
			os.mkdir(self.object_name)
		if not os.path.exists(self.object_name+'/klipped_cubes_Nfakes'):
			os.mkdir(self.object_name+'/klipped_cubes_Nfakes')

		# Determinng Number of KLIP Runs That Will Be Run
		number_of_klip = len(self.trials)

		self.write_to_log_and_print('####### BEGINNING KLIP ON DATA WITHOUT FAKES #######\n####### Number of KLIP Runs '
									'To Complete: {0} #######\n'.format(number_of_klip))

		for i, trial in enumerate(self.trials): # i only used for measuring progress
			# Get Number of KLIP Runs Conducted Already
			klip_runs = i

			# Update Every 5
			if (klip_runs + 1) % 5 == 0:
				self.write_to_log_and_print("####### {0}/{1} KLIP Runs Complete ({2}%) #######\n".
											format(klip_runs + 1, number_of_klip, round(float(klip_runs + 1) /
																					 float(number_of_klip),3)*100))

			if run_on_nofakes:
				# Running KLIP on Data Without Fakes
				with log_file_output(self.object_name):
					klip_dataset(self.dataset, outputdir=self.object_name+'/klipped_cubes_Nfakes',
								 fileprefix=self.object_name+ '_withoutfakes_' + trial.klip_parameters,
								 annuli=trial.annuli, subsections=trial.subsections, movement=trial.movement,
								 numbasis=trial.numbasis, spectrum=trial.spectrum, verbose=True,
								 corr_smooth=trial.corr_smooth, highpass=trial.highpass, mode=self.mode,
								 numthreads=65)

			# Update Every 5 or When Completely Done
			if i + 1 == len(self.trials):
				self.write_to_log_and_print("\n### DONE WITH KLIP ON DATA WITHOUT FAKES###")
			elif (klip_runs + 2) % 10 == 0:
				self.write_to_log_and_print("####### {0}/{1} KLIP Runs Complete ({2}%) #######".
											format(klip_runs + 2, number_of_klip, round(float(klip_runs + 2) /
																					 float(number_of_klip),3)*100))


	def run_KLIP_on_data_with_fakes(self):
		# Making Sure Output Directories Exist
		if not os.path.exists(self.object_name):
			os.mkdir(self.object_name)
		if not os.path.exists(self.object_name+'/klipped_cubes_Wfakes'):
			os.mkdir(self.object_name+'/klipped_cubes_Wfakes')

		number_of_klip = len(self.trials)

		self.write_to_log_and_print('####### BEGINNING KLIP ON DATA WITH FAKES #######\n####### Number of KLIP Runs To '
									'Complete: {0} #######\n'.format(number_of_klip))

		for i, trial in enumerate(self.trials): # i only used for measuring progress
			# Get Number of KLIP Runs Conducted Already
			klip_runs = i

			# Running KLIP
			with log_file_output(self.object_name):
				klip_dataset(self.dataset, outputdir=self.object_name + '/klipped_cubes_Wfakes',
							 fileprefix=self.object_name + '_withfakes_' + trial.klip_parameters,
							 annuli=trial.annuli, subsections=trial.subsections, movement=trial.movement,
							 numbasis=trial.numbasis, spectrum=trial.spectrum, verbose=True,
							 corr_smooth=trial.corr_smooth, highpass=trial.highpass, mode=self.mode,
							 numthreads=65)
				trial.output_wcs = self.dataset.output_wcs[0] # need this for calibrating contrast curve

			# Update Every 5 or When Completely Done
			if i + 1 == len(self.trials):
				self.write_to_log_and_print("\n### DONE WITH KLIP ON DATA WITH FAKES ###")
			elif (klip_runs + 1) % 5 == 0:
				self.write_to_log_and_print("####### {0}/{1} KLIP Runs Complete ({2}%) #######".
											format(klip_runs + 2, number_of_klip, round(float(klip_runs + 2) /
																					 float(number_of_klip),3)*100))


	def contrast_and_detection(self, detect_planets=True, datasetwithfakes=[True]):
		if True in datasetwithfakes or False in datasetwithfakes or detect_planets == True:
			# Making Sure Directories Exist
			if not os.path.exists(self.object_name + '/detections') and detect_planets:
				os.mkdir(self.object_name + '/detections')
			if not os.path.exists(self.object_name + '/calibrated_contrast') and True in datasetwithfakes:
				os.mkdir(self.object_name + '/calibrated_contrast')
			if not os.path.exists(self.object_name + '/uncalibrated_contrast') and False in datasetwithfakes:
				os.mkdir(self.object_name + '/uncalibrated_contrast')

			self.write_to_log_and_print("\n############## BEGINNING CONTRAST AND DETECTION FOR {0} "
										 "##############".format(self.object_name))

			for i, trial in enumerate(self.trials): # i only used for progress updates
				# if calib=True and fake planets injected, contrast will be calibrated w/ respect to KLIP subtraction
				for calib in datasetwithfakes:
					trial.get_contrast(calib)
				if detect_planets:
					# if withfakes=True, detections will use KLIP output w/ fakes, else will use KLIP output w/o fakes
					for withfakes in datasetwithfakes:
						trial.detect_planets(datasetwithfakes=withfakes)

				# Update Every 10 or When Completely Done
				if i + 1 == len(self.trials):
					self.write_to_log_and_print('\n############## DONE WITH CONTRAST AND DETECTION FOR {0} '
												 '##############'.format(self.object_name))
				elif ((i+1) * 2) % 10 == 0:
					self.write_to_log_and_print("\n####### Detection and contrast complete for {0}/{1} Trials ({2}%) "
												 "#######".format((i+1) * 2, len(self.trials) * 2,
																  round(float(i+1) / float(len(self.trials)),3)*100))
		else:
			self.write_to_log_and_print("\ncontrast_and_detection function was called, but no contrast measurements or "
				  "planet detections were conducted. Check arguments.")