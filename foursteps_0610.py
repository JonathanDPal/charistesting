from astropy.io import fits
import numpy as np
from glob import glob
from pyklip.instruments.CHARIS import CHARISData
from copy import copy


def FWHMIOWA_calculator(speccubefile, D=8, lenslet_scale=0.0162,
						field_radius=1.035):
	"""
	Finds FWHM, IWA, and OWA for a opened CHARIS data cube.
	---
	Args:
		speccubefile: Opened fits file of data cube.
		D: The aperature of the telescope, in meters. Default: 8.
		lenslet_scale: Default: 0.0162 arcsec/pixel.
		field_radius: Should be entered in arcseconds. Default: 1.035
	Returns:
		FWHM, IWA, and OWA of dataset
	"""
	wavelengths = {'j': 1200e-9, 'h': 1550e-9, 'k': 2346e-9,
				   'broadband': 1550e-9}
	wavelength = wavelengths[str.lower(speccubefile[1].header['FILTNAME'])]
	FWHM = 2 * 1.22 * wavelength / D * 206265 / lenslet_scale
	IWA = 5
	OWA = (field_radius / lenslet_scale) - FWHM

	return FWHM, IWA, OWA


def calibrate_ss_contrast(speccubefile):
	"""
	Provides ratio to calibrate flux to be relative to the star.
	---
	Args:
		speccubefile: The KLIPped datacube to get info from and for

	Returns:
		wln_um: wavelengths in micrometers
		spot_to_star: spot/star flux ratio for each wavelength
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


def step_one(fileset, fake_fluxes, fake_seps, output_filepath, fake_fwhm=3.5,
			 fake_PAs=[0,90,180,270]):
	"""
	Injects fake planets into CHARIS data and saves the data as a FITS file.
	---
	Args:
		fileset (str): Should identify where to find CHARIS data files;
					   eg. 'HD1160_cubes/*.fits'
		fake_fluxes (list of integers): The fluxes to give to injected fake
										planets.
		fake_seps (list of integers): The seperations at which to inject
									  fake planets.
		fake_fwhm (float): The fwhm for fake planet injection. Default: 3.5.
		fake_PAs (list of integers): The position angles at which to inject
									 fake planets. Default: [0,90,180,270].
		output_filepath (str): The filepath to save data with fakes injected
							   to. Should end in '.fits'
	
	Returns:
		None. 

		Will save output of injecting fake planets as an FITS file which,
		when opened, will be an HDUList object containing an empty PrimaryHDU
		in extension 0 and an ImageHDU with the data with fake planets
		injected in extension 1.
	"""
	from pyklip.fakes import inject_planet

	# Load (Actual) CHARIS Data
	filelist = glob(fileset)
	dataset = CHARISData(filelist, guess_spot_locs=None)

	# Getting Values
	with fits.open(filelist[0]) as hdu:
		_, spot_to_star = calibrate_ss_contrast(hdu)

	# Inject Fake Planets
	for fake_flux, sep in zip(fake_fluxes, fake_seps):
		flux_to_inject = fake_flux * spot_to_star
		for pa in fake_PAs:
			inject_planet(dataset.input, dataset.centers, flux_to_inject,
					    dataset.wcs, sep, pa, fwhm=fake_fwhm)

	# Save Data With Fakes Injected
	datahdu = fits.ImageHDU(dataset.input)
	datahdu.writeto(output_filepath)


def step_two(object_name, fileset, step_one_output_filepath,
			 without_fakes_output_directory, with_fakes_output_directory,
			 annuli, subsections, movement, numbasis, fake_fluxes,
			 fake_seps, fake_PAs, mask_xy=None, mode='ADI+SDI',
			 KL_modes_contrast_curve=20):
	"""
	Runs KLIP on data and creates two contrast curves and two csv files with
	the data plotted on the contrast curves. The first contrast curve and
	associated csv file will be calibrated with respect to star flux,
	but not calibrated with respect to KLIP attenuation. The second will be
	calibrated with respect to both.
	---
	Args:
		object_name (str): Identifier for data being used.
		fileset (str): Should identify where to find CHARIS data files;
					   eg. 'HD1160_cubes/*.fits'
		step_one_output_filepath (str): Filepath where output of step one
										was saved.
		without_fakes_output_directory (str): Directory to save KLIP output
											  of data without fakes to.
		with_fakes_output_directory (str): Directory to save KLIP output of
										   data with fakes to.
		annuli (int): KLIP parameter. See PyKLIP documentation.
		subsections (int): KLIP parameter.
		movement (int): KLIP parameter.
		numbasis (int or list of integers): KLIP parameter.
		fake_fluxes (list): Values used for variable in step one.
		fake_seps (list): Values used for variable in step one.
		fake_PAs (list): Values used for variable in step one.
		mask_xy (list): If a planet/companion needs to be masked, enter list
						of form [X-coor, Y-coor] giving location. Default: None.
		mode (str): KLIP parameter. 'ADI+SDI' by default.
		KL_modes_contrast_curve (int): Identifies which KLIP output to use
									   for contrast curves. 20 by default.

	Returns:
		None.
 
		Function will save KLIPped datasets as FITS files, contrast curves
		as PNG images, and contrast vs. seperation data as CSV files in
		respective directories.
	"""
	### Imports ###
	from pyklip.parallelized import klip_dataset
	from pyklip.klip import meas_contrast

	### Helper Functions For Step Two ###

	def ss_calibration(contains_fakes):
		"""
		Only intended for use within step_two function; won't work outside
		of it. Calibrates data cube with respect to star flux and provides a
		bunch of other information needed for contrast measurements.
		---
		Arg:
			contains_fakes (bool): Set to True if looking to calibrate data
								   with the fakes and to False if looking to
								   calibrate data without the fakes.

		Returns:
			wavelength in micrometers
			datacube calibrated with respect to star flux
			dataset center as a list of form [X-coor, Y-coor]
			dataset inner working angle
			dataset outer working angle
			dataset full width at half maximum
			dataset output WCS
		"""
		from astropy.wcs import WCS


		if contains_fakes:
			directory = with_fakes_output_directory
			object_name_end = '_withfakes-KL'
		elif not contains_fakes:
			directory = without_fakes_output_directory
			object_name_end = '_withoutfakes-KL'

		with fits.open(directory+'/'+object_name+object_name_end
				+str(KL_modes_contrast_curve)+'-speccube.fits') as hdulist:
			wln_um, spot_to_star = calibrate_ss_contrast(hdulist)
			calib_cube = copy(hdulist[1].data) * spot_to_star[:, np.newaxis,
												 np.newaxis]
			dataset_center = [hdulist[1].header['PSFCENTX'],
							  hdulist[1].header['PSFCENTY']]
			dataset_fwhm,dataset_iwa,dataset_owa = FWHMIOWA_calculator(hdulist)
			output_wcs = WCS(header=hdulist[1].header, naxis=[1,2])

		return wln_um, calib_cube, dataset_center, dataset_iwa, dataset_owa, \
			   dataset_fwhm, output_wcs


	def generate_contrast_curve(graph_output_filepath, data_output_filepath,
								wln_um, calib_cube, dataset_center,
								dataset_iwa, dataset_owa, dataset_fwhm,
								output_wcs, wavelength_index=10,
								mask_xy=mask_xy, contains_fakes=False,
								injected_fluxes=None, injected_seps=None,
								injected_PAs=None):
		"""
		Measures contrast in an image and generates contrast curve.
		---
		Args:
			graph_output_filepath (str): Filepath for contrast curve image.
										 Should end in '.png'
			data_output_filepath (str): Filepath for contrast curve data.
										Should end in '.csv'
			wln_um (list): List of wavelengths in micrometers.
			calib_cube (np array): Data cube which has been calibrated with
								   regard to star/spot flux.
			dataset_center (list)
			dataset_iwa (float or int)
			dataset_owa (float)
			dataset_fwhm (float)
			output_wcs (instance of astropy.wcs WCS class)
			wavelength_index (int): Index of wavelength to use (subsets
									calibrated cube along wavelength axis).
									Default: 10
			mask_xy (list): Reading in argument from step_two func.
			contains_fakes (bool): Set to true if fakes present. Default: False
			injected_fluxes (list): If fakes present, provide list of fake
									planets' fluxes. Default: None
			injected_seps (list): If fakes present, provide list of fake
								  planets' seperations. Default: None
			injected_PAs (list): If fakes present, provide list of fake
								 planets' position angles. Default: None


		Returns:


		"""
		import pylab as P
		from csv import writer

		frame = calib_cube[wavelength_index]

		if mask_xy is not None:
			x_pos = mask_xy[0]
			y_pos = mask_xy[1]

			ydat, xdat = np.indices(frame.shape)
			distance_from_planet = np.sqrt((xdat-x_pos)**2+(ydat-y_pos)**2)
			frame[np.where(distance_from_planet <= 2*dataset_fwhm)] = np.nan

		contrast_seps, contrast = meas_contrast(frame, dataset_iwa,
												dataset_owa, dataset_fwhm,
												center=dataset_center,
												low_pass_filter=True)

		# Calibrating For KLIP Subtraction If Fakes Present
		if contains_fakes:
			from pyklip.fakes import retrieve_planet_flux
			retrieved_fluxes = []
			for sep in injected_seps:
				fake_planet_fluxes = []
				for pa in injected_PAs:
					fake_flux = retrieve_planet_flux(frame, dataset_center,
													 output_wcs, sep, pa,
													 searchrad=7)
					fake_planet_fluxes.append(fake_flux)
				retrieved_fluxes.append(np.mean(fake_planet_fluxes))

			algo_throughput = \
				np.array(retrieved_fluxes) / np.array(injected_fluxes)

			correct_contrast = np.copy(contrast)
			for i, sep in enumerate(contrast_seps):
				closest_throughput_index = np.argmin(np.abs(sep-injected_seps))
				correct_contrast[i] /= algo_throughput[closest_throughput_index]

		# Plotting & Saving Curve
		if contains_fakes:
			P.semilogy(contrast_seps, correct_contrast)
		elif not contains_fakes:
			P.semilogy(contrast_seps, contrast)
		P.xlabel('Seperation (pixels)')
		if contains_fakes:
			P.ylabel(r'$5\sigma$ Contrast (Calibrated)')
		elif not contains_fakes:
			P.ylabel(r'$5\sigma$ Contrast (Not Calibrated)')
		wavelength = round(wln_um[wavelength_index], 2)
		if contains_fakes:
			P.title(r'$5\sigma$ Contrast at {0} $\mu$m (Calibrated)'.format(
				wavelength))
		elif not contains_fakes:
			P.title(r'$5\sigma$ Contrast at {0} $\mu$m (Not '
					r'Calibrated)'.format(wavelength))
		P.savefig(graph_output_filepath)

		# Saving Data
		if contains_fakes:
			with open(data_output_filepath, 'w+') as csvfile:
				csvwriter = writer(csvfile, delimiter=',')
				csvwriter.writerows([['Sep (Pixels)', 'Contrast']])
				csvwriter.writerows([contrast_seps, correct_contrast])
		elif not contains_fakes:
			with open(data_output_filepath, 'w+') as csvfile:
				csvwriter = writer(csvfile, delimiter=',')
				csvwriter.writerows([['Sep (Pixels)', 'Contrast']])
				csvwriter.writerows([contrast_seps, contrast])

	##############################################
	########## Actual Start of Step Two ##########
	##############################################

	# Loading Data
	filelist = glob(fileset)
	data_without_fakes = CHARISData(filelist, guess_spot_locs=None)
	data_with_fakes = copy(data_without_fakes)
	fakes_data = fits.getdata(step_one_output_filepath, 1)
	data_with_fakes.input = fakes_data # now it has fake planets from step one

	# Running KLIP on Data Without Fakes
	klip_dataset(data_without_fakes,
				 outputdir=without_fakes_output_directory,
				 fileprefix=object_name+'_withoutfakes', annuli=annuli,
				 subsections=subsections, movement=movement,
				 numbasis=numbasis, mode=mode)

	# Measuring Contrast, Generating Contrast Curve on Data Without Fakes &
	# Saving the Contrast Data to a CSV File
	no_fakes_wln_um, no_fakes_calib_cube, no_fakes_dataset_center, \
	no_fakes_dataset_iwa, no_fakes_dataset_owa, no_fakes_dataset_fwhm, \
	no_fakes_output_wcs = ss_calibration(False)

	filepath_uncalibrated_curve = '' # where contrast curve will be saved
	filepath_uncal_contrast_data = '' # where contrast data will be saved

	generate_contrast_curve(graph_output_filepath=filepath_uncalibrated_curve,
							data_output_filepath=filepath_uncal_contrast_data,
							wln_um=no_fakes_wln_um,
							calib_cube=no_fakes_calib_cube,
							dataset_center=no_fakes_dataset_center,
							dataset_iwa=no_fakes_dataset_iwa,
							dataset_owa=no_fakes_dataset_owa,
							dataset_fwhm=no_fakes_dataset_fwhm,
							output_wcs=no_fakes_output_wcs)

	# Running KLIP on Data With Fakes
	klip_dataset(data_with_fakes, outputdir=with_fakes_output_directory,
				 fileprefix=object_name+'_withfakes', annuli=annuli,
				 subsections=subsections, movement=movement,
				 numbasis=numbasis, mode=mode)

	# Measuring Contrast, Generating Contrast Curve on Data With Fakes &
	# Saving the Contrast Data to a CSV File
	with_fakes_wln_um, with_fakes_calib_cube, with_fakes_dataset_center, \
	with_fakes_dataset_iwa, with_fakes_dataset_owa, with_fakes_dataset_fwhm, \
	with_fakes_output_wcs = ss_calibration(True)

	filepath_calibrated_curve = '' # where contrast curve will be saved
	filepath_cal_contrast_data = '' # where contrast data will be saved

	generate_contrast_curve(graph_output_filepath=filepath_calibrated_curve,
							data_output_filepath=filepath_cal_contrast_data,
							wln_um=with_fakes_wln_um,
							calib_cube=with_fakes_calib_cube,
							dataset_center=with_fakes_dataset_center,
							dataset_iwa=with_fakes_dataset_iwa,
							dataset_owa=with_fakes_dataset_owa,
							dataset_fwhm=with_fakes_dataset_fwhm,
							output_wcs=with_fakes_output_wcs,
							contains_fakes=True, injected_fluxes=fake_fluxes,
							injected_seps=fake_seps, injected_PAs=fake_PAs)


def step_three(filename, output_prefix, SNR_threshold=3):
	"""
	Looks at a KLIPped dataset and indicates potential planets.
	---
	Args:
		filename: Path to the KLIPped dataset to be used.
		output_prefix: Beginning of filepath for output. Output will have the
					   form "(output_prefix)_SNR-(SNR_threshold).csv".
		SNR_threshold: Default: 3
	"""
	from pyklip.kpp.utils.mathfunc import gauss2d
	from pyklip.kpp.metric.crossCorr import calculate_cc
	from pyklip.kpp.stat.statperPix_utils import \
		get_image_stat_map_perPixMasking
	from pyklip.kpp.detection.detection import point_source_detection
	from csv import writer

	with fits.open(filename) as hdulist:
		image = hdulist[1].data
		center = [hdulist[0].header['PSFCENTX'], hdulist[0].header['PSFCENTY']]

	x_grid, y_grid = np.meshgrid(np.arange(-10,10), np.arange(-10,10))
	kernel_gauss = gauss2d(x_grid, y_grid)

	image_cc = calculate_cc(image, kernel_gauss,
							spectrum=np.ones(len(image[0])), nans2zero=True)

	SNR_map = get_image_stat_map_perPixMasking(image_cc, centroid=center,
											   mask_radius=5, Dr=2, type='SNR')

	detection_threshold = SNR_threshold
	pix2as = 1
	mask_radius = 15
	maskout_edge = 10

	candidates_table = point_source_detection(SNR_map, center,
											  detection_threshold,
											  pix2as=pix2as,
											  mask_radius=mask_radius,
											  maskout_edge=maskout_edge,
											  IWA=None, OWA=None)

	# Write Out Detections To CSV File
	with open(output_prefix+'_SNR-'+str(SNR_threshold)+'.csv','w+') as csvfile:
		csvwriter = writer(csvfile, delimiter=',')
		csvwriter.writerows([['Index', 'SNR Value', 'PA', 'Sep (pix)',
							  'Sep (as)', 'x', 'y', 'row', 'col']])
		csvwriter.writerows(candidates_table)


def step_four(detections_finders, contrasts):
	"""
	Aggregates information and provides summary plots.
	---
	Args:
		detections_finders (string or list of strings): String(s) will be
														passed into glob in
														order to find the csv
														files with detections.
		contrasts (string or list of strings): Strings(s) will be passed into
											   glob in order to find the csv
											   files with the contrast vs.
											   seperation data.
	"""
	import pandas as pd

	detections_list = []
	for detection_finder in detections_finders:
		detections_list.append(glob(detection_finder))
	detections = dict()
	for filename in detections_list:
		if filename[-6] == '-':
			snr_value = filename[-5]
		elif filename[-7] == '-':
			snr_value = filename[-6] + filename[-5]
		else:
			print("Filename {0} is not correctly formatted and will not be "
				  "included in the compiled data")
			continue
		detections[snr_value] = pd.read_csv(filename)