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


def calibrate_ss_contrast(speccubefile):
	"""
	Provides ratios to calibrate flux to be relative to the star at all wavelengths. Thanks to
	Dr. Tobin for this. (https://docs.google.com/document/d/1S1Oo9QweKwnOfmv6fu28bb75lYeQXGzn/edit)
	"""
	cube = deepcopy(speccubefile[1].data)

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


def pasep_to_xy(PAs, seps):
	"""
	Takes lists of position angles and seperations and yields a numpy array with x-y coordinates for each combo.
	"""
	radians = np.array(PAs) / 180 * np.pi
	xy = []
	for sep in seps:
		for rad in radians:
			x = -np.sin(rad) * sep
			y = np.cos(rad) * sep
			xy.append(np.array([x,y]))

	return np.array(xy)


# TestDataset Will Have a List of Trials Associated With It (one for each group of KLIP Parameters)
class Trial:
	"""
	NOTE: The user will almost certainly not interact with this class directly, rather they will interact with an
	instance of TestDataset and that instance of TestDataset will interact with instances of this class.
	---
	Stores a particular set of KLIP parameters and then is able to run contrast measurement or planet detection code
	from KLIP for the KLIP output with that particular set of parameters.
	"""
	def __init__(self, annuli, subsections, movement, numbasis, spectrum, mask_xy, fake_PAs, fake_fluxes,
				 object_name, fake_fwhm, fake_seps, corr_smooth, highpass, length):
		self.annuli = annuli
		self.subsections = subsections
		self.movement = movement
		self.numbasis = numbasis
		self.spectrum = spectrum
		self.mask_xy = mask_xy
		self.fake_PAs = np.array(fake_PAs)
		self.fake_fluxes = fake_fluxes
		self.object_name = object_name
		self.fake_fwhm = fake_fwhm
		self.fake_seps = np.array(fake_seps)
		self.corr_smooth = corr_smooth

		# Switching Highpass From Fourier to Pixels If Necessary
		if isinstance(highpass, (int, float)) and not isinstance(highpass, bool):
			highpass = float(highpass)
			self.highpass = length / (highpass * 2 * np.sqrt(2 * np.log(2)))
		else:
			self.highpass = highpass

		# String Identifying Parameters Used
		self.klip_parameters = str(annuli)+'Annuli_'+str(subsections)+'Subsections_'+str(movement)+'Movement_'+str(
								spectrum)+'Spectrum_'+str(corr_smooth)+'Smooth_'+str(highpass)+'Highpass_'

		# Filepaths to KLIPped Datacubes
		self.filepaths_Wfakes = [self.object_name + '/klipped_cubes_Wfakes/' + self.object_name + '_withfakes_' +
								 self.klip_parameters +'-KL{0}-speccube.fits'.format(nb) for nb in self.numbasis]
		self.filepaths_Nfakes = [self.object_name + '/klipped_cubes_Nfakes/' + self.object_name + '_withoutfakes_' +
								 self.klip_parameters +'-KL{0}-speccube.fits'.format(nb) for nb in self.numbasis]

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
		Measures contrast, then saves a contrast curve as a PNG and the contrast data as a CSV.
		---
		Args:
			contains_fakes (bool): Set to True if fakes present.
			wavelength_index (int): Index of wavelength to use (subsets calibrated cube along wavelength axis).
		"""
		if contains_fakes:
			filepaths = self.filepaths_Wfakes
		else:
			filepaths = self.filepaths_Nfakes

		for i, filepath in enumerate(filepaths):
			with fits.open(filepath) as hdulist:
				wln_um, spot_to_star = calibrate_ss_contrast(hdulist)
				calib_cube = deepcopy(hdulist[1].data) * spot_to_star[:, np.newaxis, np.newaxis]
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

			contrast_seps, contrast = meas_contrast(frame, dataset_iwa, dataset_owa, dataset_fwhm,
													center=dataset_center, low_pass_filter=True)

			# Calibrating For KLIP Subtraction If Fakes Present
			if contains_fakes:
				retrieved_fluxes = []
				for sep in self.fake_seps:
					fake_planet_fluxes = []
					for pa in self.fake_PAs:
						fake_flux = retrieve_planet_flux(frame, dataset_center, output_wcs, sep, pa, searchrad=7)
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
			if contains_fakes:
				if not os.path.exists(self.object_name + '/calibrated_contrast'):
					os.mkdir(self.object_name + '/calibrated_contrast')
			else:
				if not os.path.exists(self.object_name + '/uncalibrated_contrast'):
					os.mkdir(self.object_name + '/uncalibrated_contrast')

			# Saving Data to Object and to CSV File
			if contains_fakes:
				data_output_filepath = self.object_name + '/calibrated_contrast/{0}_KL{1}_contrast.csv'.format(
										   self.klip_parameters, self.numbasis[i])
				df = pd.DataFrame()
				df['Seperation'] = contrast_seps
				df['Calibrated Contrast'] = contrast
				wavelength = wln_um[wavelength_index]
				title = 'Calibrated Contrast at {0}um'.format(wavelength)
				df.plot(x='Seperation', y='Calibrated Contrast', legend=False, title=title)
				plt.ylabel('Calibrated Contrast')
				plt.xlabel('Seperation')
				plt.semilogy()
				plt.savefig(data_output_filepath[0:-4] + '.png')
				df.to_csv(data_output_filepath)
			else:
				data_output_filepath = self.object_name + '/uncalibrated_contrast/{0}_KL{1}_contrast.csv'.format(
										   self.klip_parameters, self.numbasis[i])
				df = pd.DataFrame()
				df['Seperation'] = contrast_seps
				df['Uncalibrated Contrast'] = contrast
				wavelength = round(wln_um[wavelength_index], 2)
				title = 'Uncalibrated Contrast at {0}um ({1})'.format(wavelength, self.object_name)
				df.plot(x='Seperation', y='Uncalibrated Contrast', legend=False, title=title)
				plt.ylabel('Uncalibrated Contrast')
				plt.xlabel('Seperation')
				plt.semilogy()
				plt.savefig(data_output_filepath[0:-4] + '.png')
				df.to_csv(data_output_filepath)


	def detect_planets(self, SNR_threshold=3):
		"""
		Looks at a KLIPped dataset with fakes and indicates potential planets.
		---
		Args:
			SNR_threshold: Default: 3. Set this to the lowest value to be looked at; when
			generating ROC curves, the code will automatically look at subsets of the detections
			at different SNR thresholds, but it can only look at SNR thresholds greater than or
			equal to the one specifie here.
		"""
		for i, filepath in enumerate(self.filepaths_Wfakes):
			with fits.open(filepath) as hdulist:
				image = hdulist[1].data
				center = [hdulist[1].header['PSFCENTX'], hdulist[1].header['PSFCENTY']]

			x_grid, y_grid = np.meshgrid(np.arange(-10,10), np.arange(-10,10))
			kernel_gauss = gauss2d(x_grid, y_grid)

			# flat spectrum given so that it collapses it into one image, instead of giving
			# seperate images for each wavelength
			image_cc = calculate_cc(image, kernel_gauss, spectrum=np.ones(len(image[0])), nans2zero=True)

			SNR_map = get_image_stat_map_perPixMasking(image_cc, centroid=center, mask_radius=5, Dr=2, type='SNR')

			candidates_table = point_source_detection(SNR_map, center, SNR_threshold, pix2as=1, mask_radius=15,
													  maskout_edge=10, IWA=None, OWA=None)

			candidates = pd.DataFrame(candidates_table, columns=['Index', 'SNR Value', 'PA', 'Sep (pix)',
																 'Sep (as)', 'x', 'y', 'row', 'col'])

			injected = []
			fakelocs = pasep_to_xy(self.fake_PAs, self.fake_seps)
			for _, row in candidates.iterrows():
				distances = [np.sqrt((row['x'] - fl[0]) ** 2 + (row['y'] - fl[1]) ** 2) for fl in fakelocs]
				if np.min(distances) > self.fake_fwhm:
						injected.append(False)
				else:
						injected.append(True)
			candidates['Injected'] = injected
			candidates.to_csv('{0}{1}.csv'.format(self.filepath_detections_prefixes[i], str(SNR_threshold)))


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

		print("##################### STARTING WORK ON {0} #####################".format(self.object_name))

		# Creating CHARISData Object With UnKLIPped Data
		self.fileset = glob(fileset)
		self.dataset_no_fakes = CHARISData(self.fileset, guess_spot_locs=[[91,72], [129, 91], [71, 111], [109, 130]],
										   guess_center_loc=[100,100])
		print("####### DONE BUILDING CHARISData OBJECT FOR {0} ########".format(self.object_name))
		self.dataset_with_fakes = deepcopy(self.dataset_no_fakes)
		self.length = self.dataset_no_fakes.input.shape[1]

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
									self.trials.append(Trial(annuli=ani, subsections=subsec,movement=mov,
															 numbasis=numbasis, spectrum=spec,mask_xy=mask_xy,
															 fake_PAs=fake_PAs,fake_fluxes=fake_fluxes,
															 object_name=object_name,fake_fwhm=fake_fwhm,
															 fake_seps=fake_seps, corr_smooth=cs,
															 highpass=hp, length=self.length))
		self.mode = mode
		print("############## DONE BUILDING TRIALS FOR {0} ##############".format(self.object_name))

	def inject_fakes(self):
		# Getting Values
		with fits.open(self.fileset[0]) as hdu:
			spot_to_star = calibrate_ss_contrast(hdu)[1]

		# Inject Fake Planets
		for fake_flux, sep in zip(self.fake_fluxes, self.fake_seps):
			flux_to_inject = fake_flux / spot_to_star # UNcalibrating it, NOT calibrating
			for pa in self.fake_PAs:
				inject_planet(frames=self.dataset_with_fakes.input, centers=self.dataset_with_fakes.centers,
							  inputflux=flux_to_inject, astr_hdrs=self.dataset_with_fakes.wcs, radius=sep, pa=pa,
							  fwhm=self.fake_fwhm)

		print("############## DONE INJECTING FAKES FOR {0} ##############".format(self.object_name))


	def run_KLIP(self, run_on_fakes=True, run_on_nofakes=True):
		if run_on_fakes or run_on_nofakes: # making sure at least set of KLIP runs will happen
			# Making Sure Output Directories Exist
			if not os.path.exists(self.object_name):
				os.mkdir(self.object_name)
			if run_on_fakes:
				if not os.path.exists(self.object_name+'/klipped_cubes_Wfakes'):
					os.mkdir(self.object_name+'/klipped_cubes_Wfakes')
			if run_on_nofakes:
				if not os.path.exists(self.object_name+'/klipped_cubes_Nfakes'):
					os.mkdir(self.object_name+'/klipped_cubes_Nfakes')

			if run_on_fakes and run_on_nofakes:
				running_on_both = True
			else:
				running_on_both = False

			# Determining Number of KLIP Runs That Will Be Conducted
			if running_on_both:
				number_of_klip = len(self.trials) * 2
			else:
				number_of_klip = len(self.trials)

			print("############## BEGINNING KLIP FOR {0} ##############".format(self.object_name))
			print("####### Total KLIP Runs to Complete: {0} #######".format(number_of_klip))

			for i, trial in enumerate(self.trials):
				if running_on_both:
					trials = i * 2
				else:
					trials = i

				if run_on_fakes:
					# Running KLIP on Data With Fakes
					with log_file_output():
						klip_dataset(self.dataset_with_fakes, outputdir=self.object_name+'/klipped_cubes_Wfakes',
									 fileprefix=self.object_name + '_withfakes_' + trial.klip_parameters,
									 annuli=trial.annuli, subsections=trial.subsections, movement=trial.movement,
									 numbasis=trial.numbasis, spectrum=trial.spectrum, verbose=False,
									 corr_smooth=trial.corr_smooth, highpass=trial.highpass, mode=self.mode)

				# Update Every 5
				if (trials + 1) % 5 == 0:
					print("####### {0}/{1} KLIP Runs Complete ({2}%) #######".format(trials + 1, number_of_klip,
																					 round(float(trials + 1) / float(
																						 number_of_klip), 3) * 100))

				if run_on_nofakes:
					# Running KLIP on Data Without Fakes
					with log_file_output():
						klip_dataset(self.dataset_no_fakes, outputdir=self.object_name+'/klipped_cubes_Nfakes',
									 fileprefix=self.object_name+ '_withoutfakes_' + trial.klip_parameters,
									 annuli=trial.annuli, subsections=trial.subsections, movement=trial.movement,
									 numbasis=trial.numbasis, spectrum=trial.spectrum, verbose=False,
									 corr_smooth=trial.corr_smooth, highpass=trial.highpass, mode=self.mode)

				# Update Every 5 or When Completely Done
				if i + 1 == len(self.trials):
					print("############## DONE WITH KLIP FOR {0} ##############".format(self.object_name))
				elif (trials + 2) % 10 == 0:
					print("####### {0}/{1} KLIP Runs Complete ({2}%) #######".format(trials + 2, number_of_klip,
																					 round(float(trials + 2) / float(
																						 number_of_klip), 3) * 100))
		else:
			print("run_KLIP function called, but no KLIP runs conducted. Check arguments.")


	def contrast_and_detection(self, calibrate=[True, False], detect_planets=True):
		if True in calibrate or False in calibrate or detect_planets == True: # checking arguments
			print("############## BEGINNING CONTRAST AND DETECTION FOR {0} ##############".format(self.object_name))
			for i, trial in enumerate(self.trials):
				for calib in calibrate:
					trial.get_contrast(calib)
				if detect_planets:
					trial.detect_planets()
				if i + 1 == len(self.trials):
					print('############## DONE WITH CONTRAST AND DETECTION FOR {0} ##############'.format(
						self.object_name))
				elif ((i+1) * 2) % 10 == 0:
					print("####### Detection and contrast complete for {0}/{1} Trials ({2}%) #######".format((i+1) *
							2, len(self.trials) * 2, round(float(i+1) / float(len(self.trials)), 3) * 100))
		else:
			print("contrast_and_detection function was called, but no contrast measurements or "
				  "planet detections were conducted. Check arguments.")


	@contextmanager
	def log_file_output(self):
		with open('{0}/log.txt'.format(self.object_name), 'w') as log_file:
			old_stdout = sys.stdout
			sys.stdout = log_file
			try:
				yield
			finally:
				sys.stdout = old_stdout


	def write_to_log(self, words):
		with open('{0}/log.txt'.format(self.object_name), 'w') as log_file:
			log_file.write(words)
