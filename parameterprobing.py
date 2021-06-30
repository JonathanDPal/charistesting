from parameter_test_infrastructure import *
import warnings
from time import time
from numpy import floor
from astropy.io import fits
from glob import glob

#####################################################
#################### USER INPUTS ####################
#####################################################

# see https://docs.google.com/document/d/1yX0l96IZs1IxxKCRmriVSAQM3KFGF9U1-FnpJXhcLXo/edit?usp=sharing for help

# General Set-Up
fileset0 = 'HD1160_cubes_V2/*.fits'
mask0 = None
object_name0 = 'HD1160_new_dn_per_contrast'

# Setting Up For KLIP
# annuli = [4, 6, 8, 10, 12]
annuli = [9]
# subsections = [2, 4, 6]
subsections = [4]
movement = [1]
spectrum = [None]
numbasis = [20]
# numbasis = [10, 20, 30, 40, 50, 60]
# corr_smooth = [0, 1, 2]
corr_smooth = [1]
highpass = [True]
# highpass = [False, 5.0, True, 15.0] # True yields default (10.0)
mode = 'ADI+SDI'

# Setting Up For Fake Planets
fake_fluxes = [1e-2, 1e-3, 1e-4]
fake_seps = [15, 30, 45, 60]
fake_PAs=[0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]

# Specifying Particular Things To (Not) Do
put_in_fakes = True
run_KLIP_on_dataset_with_fakes = True # if no fakes are injected, this will just be a duplicate
run_KLIP_on_dataset_without_fakes = False
get_uncalibrated_contrast = False
get_calibrated_contrast = True # won't be calibrated if no fake planets were injected
get_planet_detections = True

############################################################
#################### END OF USER INPUTS ####################
############################################################

# Synthesizing User Inputs Into a Couple Additional Booleans
get_contrast_and_detections = (get_uncalibrated_contrast or get_calibrated_contrast or get_planet_detections)
run_KLIP_on_dataset = (run_KLIP_on_dataset_with_fakes or run_KLIP_on_dataset_with_fakes)
calibrate = []
if get_calibrated_contrast:
    calibrate.append(True)
if get_uncalibrated_contrast:
    calibrate.append(False)

################## STARTING ACTUAL TESTING #################
start = time()

# KLIP yields a bunch of RuntimeWarnings that we don't need to worry about
warnings.simplefilter('ignore', category=RuntimeWarning)

with fits.open(glob(fileset0)[0]) as hdulist:
   fake_fwhm = FWHMIOWA_calculator(speccubefile=hdulist)[0]

# Create TestDataset For Each Set of Observations
td0 = TestDataset(fileset=fileset0, object_name=object_name0, mask_xy=mask0, fake_fluxes=fake_fluxes,
                  fake_seps=fake_seps, annuli=annuli, subsections=subsections, movement=movement, numbasis=numbasis,
                  corr_smooth=corr_smooth, highpass=highpass, spectrum=spectrum, mode=mode, fake_PAs=fake_PAs,
                  fake_fwhm=fake_fwhm, nofakes=run_KLIP_on_dataset_without_fakes)

# Have TestDataset Run Each Part
if put_in_fakes:
    td0.inject_fakes()
if run_KLIP_on_dataset:
    td0.run_KLIP(run_on_fakes=run_KLIP_on_dataset_with_fakes, run_on_nofakes=run_KLIP_on_dataset_without_fakes)
if get_contrast_and_detections:
    td0.contrast_and_detection(calibrate=calibrate, detect_planets=get_planet_detections)

# Print Out Time Taken
end = time()
time_elapsed = end - start
hours = int(floor(time_elapsed / 3600))
remaining_time = time_elapsed - (hours * 3600)
minutes = int(floor(remaining_time / 60))
seconds = round(remaining_time - minutes * 60)

td0.write_to_log_and_print("##################### TIME ELAPSED: {0} Hours, {1} Minutes, {2} Seconds "
                            "#####################".format(hours, minutes, seconds))


