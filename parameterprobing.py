from parameter_test_infrastructure import *
import warnings
from time import time
from numpy import floor
from astropy.io import fits
from glob import glob
import os

#####################################################
#################### USER INPUTS ####################
#####################################################

# see https://docs.google.com/document/d/1yX0l96IZs1IxxKCRmriVSAQM3KFGF9U1-FnpJXhcLXo/edit?usp=sharing for help

# General Set-Up
fileset0 = 'HD1160_cubes/*.fits' # will be passed into glob.glob() to identify extracted cubes
mask0 = [144, 80] # [X-coor, Y-coor] of object in KLIP output
object_name0 = 'HD1160' # Name of Directory Where All Outputs Will Be Placed

fileset1 = 'HR8799_cubes/*.fits'
mask1 = [[128, 152], [157, 137], [129, 70]]
object_name1 = 'HR8799'

# Setting Up Lists/Tuples For KLIP
annuli = [4, 6, 8, 10, 12] # List of Integer(s)
subsections = [2, 4, 6] # List of Integer(s)
movement = [0, 1, 2] # List or Tuple of Integer(s)
spectrum = [None] # List or Tuple of None and/or 'methane'
numbasis = [10, 20, 30, 40, 50, 60] # List or Tuple of Integer(s)
corr_smooth = [0, 1, 2] # List or Tuple of Float(s) and/or Integer(s)
highpass = [False, 5.0, True, 15.0] # List or Tuple of Float(s), Integer(s), and/or Bool(s)

# Setting Mode For KLIP
mode = 'ADI+SDI' # Exactly ONE (not a list or tuple) or the following: 'ADI', 'SDI', 'ADI+SDI'

# Setting Up For Fake Planets
fake_fluxes = [1e-3, 1e-4, 1e-5, 5e-4, 5e-5, 5e-6] # List of Float(s)
fake_seps = [20, 40, 60] # List of Integer(s) and/or Float(s)
fake_PAs=[19, 79, 139, 199, 259, 319] # List of Integer(s) and/or Float(s)

##### Specifying Which Things to Do/Not Do #####
# Most of the time, the four values below should be set to True
put_in_fakes = True
run_KLIP_on_dataset_with_fakes = True # if no fakes are injected, this will just be a dataset without fakes
get_contrast = True # won't be calibrated if no fake planets are injected
get_planet_detections_from_dataset_with_fakes = True
# Most of the time, these two values below should be set to False
run_KLIP_on_dataset_without_fakes = False
get_planet_detections_from_dataset_without_fakes = False

############################################################
#################### END OF USER INPUTS ####################
############################################################


################## CHECKING USER INPUTS #################

# Warning User if the Directory Where Stuff Will Be Outputted to Already Exists
if os.path.exists(object_name0):
    warnings.warn("WARNING: There is already a directory with the same name as the one you specified for outputs to "
                  "be written to. Outputs will overwrite previous outputs if filenames are identical")

# Making Sure This Group of Parameters Are In The Form of a List or Tuple
for param in [[annuli, 'annuli'], [subsections, 'subsections'], [movement, 'movement'], [spectrum, 'spectrum'],
              [corr_smooth, 'corr_smooth'], [highpass, 'highpass']]:
    if not isinstance(param[0], (list, tuple)):
        raise TypeError("{0} needs to be a list or a tuple. Check input.".format(param[1]))

# Checking Mode -- Common Mistake is Inputting Mode as a List/Tuple Like Other Params
if not isinstance(mode, str):
    warnings.warn("WARNING: Inputted mode is not a string. If inputted mode is a list or tuple, then script"
                  "will proceed using index 0 value as the mode. Any other values will not be used.")
    try:
        mode = mode[0]
    except Exception:
        raise TypeError("Mode needs to be a string. Check input")

# Checking Fake Planet Stuff
for param in [[fake_fluxes, 'fake_fluxes'], [fake_seps, 'fake_seps'], [fake_PAs, 'fake_PAs']]:
    if isinstance(param[0], (list, tuple)):
        if not len(fake_fluxes) % len(fake_seps) == 0:
            raise ValueError("The lengths of fake_fluxes and fake_seps must be the same.")
    else:
        if put_in_fakes:
            raise ValueError("put_in_fakes is set to true, but {0} is not a list or tuple.".format(param[1]))

##### SYNTHESIZING USER INPUTS INTO A COUPLE ADDITIONAL BOOLEANS ######
detect_planets = get_planet_detections_from_dataset_with_fakes or get_planet_detections_from_dataset_without_fakes
datasetwithfakes = []
if get_planet_detections_from_dataset_with_fakes:
    datasetwithfakes.append(True)
if get_planet_detections_from_dataset_without_fakes:
    datasetwithfakes.append(False)
get_contrast_and_detections = get_contrast or detect_planets

############################################################
################ ACTUALLY STARTING TESTING #################
############################################################
start0 = time()

# KLIP creates a bunch of RuntimeWarnings that we don't want to spam our log file
warnings.simplefilter('ignore', category=RuntimeWarning)

# FWHM Is Based on Central Wavelength of Filter Used During Obsveration.
with fits.open(glob(fileset0)[0]) as hdulist:
   fake_fwhm0 = FWHMIOWA_calculator(hdulist)[0]

# Create TestDataset For Target 0
td0 = TestDataset(fileset=fileset0, object_name=object_name0, mask_xy=mask0, fake_fluxes=fake_fluxes,
                  fake_seps=fake_seps, annuli=annuli, subsections=subsections, movement=movement, numbasis=numbasis,
                  corr_smooth=corr_smooth, highpass=highpass, spectrum=spectrum, mode=mode, fake_PAs=fake_PAs,
                  fake_fwhm=fake_fwhm0)

# Have TestDataset 0 Run Each Part
# if we want KLIP output of data without fakes, we need to run KLIP before injecting planets
if run_KLIP_on_dataset_without_fakes:
    td0.run_KLIP_on_data_without_fakes()
if put_in_fakes:
    td0.inject_fakes()
if run_KLIP_on_dataset_with_fakes:
    td0.run_KLIP_on_data_with_fakes()
if get_contrast_and_detections:
    td0.contrast_and_detection(detect_planets=detect_planets, datasetwithfakes=datasetwithfakes)

# Print Out Time Taken
end0 = time()
time_elapsed = end0 - start0
hours = int(floor(time_elapsed / 3600))
remaining_time = time_elapsed - (hours * 3600)
minutes = int(floor(remaining_time / 60))
seconds = round(remaining_time - minutes * 60)

td0.write_to_log_and_print("##################### TIME ELAPSED: {0} Hours, {1} Minutes, {2} Seconds "
                            "#####################".format(hours, minutes, seconds))

td0.dataset = None # taking HD1160 data out of memory
for t in td0.trials:
    t = None
td0 = None

#### SECOND TESTDATASET ####
start1 = time()

# FWHM Is Based on Central Wavelength of Filter Used During Obsveration.
with fits.open(glob(fileset1)[0]) as hdulist:
   fake_fwhm1 = FWHMIOWA_calculator(hdulist)[0]

# Create TestDataset For Target 1
td1 = TestDataset(fileset=fileset1, object_name=object_name1, mask_xy=mask1, fake_fluxes=fake_fluxes,
                  fake_seps=fake_seps, annuli=annuli, subsections=subsections, movement=movement, numbasis=numbasis,
                  corr_smooth=corr_smooth, highpass=highpass, spectrum=spectrum, mode=mode, fake_PAs=fake_PAs,
                  fake_fwhm=fake_fwhm1)

# Have TestDataset 1 Run Each Part
if run_KLIP_on_dataset_without_fakes:
    td1.run_KLIP_on_data_without_fakes()
if put_in_fakes:
    td1.inject_fakes()
if run_KLIP_on_dataset_with_fakes:
    td1.run_KLIP_on_data_with_fakes()
if get_contrast_and_detections:
    td1.contrast_and_detection(detect_planets=detect_planets, datasetwithfakes=datasetwithfakes)

# Print Out Time Taken
end1 = time()
time_elapsed = end1 - start1
hours = int(floor(time_elapsed / 3600))
remaining_time = time_elapsed - (hours * 3600)
minutes = int(floor(remaining_time / 60))
seconds = round(remaining_time - minutes * 60)

td1.write_to_log_and_print("##################### TIME ELAPSED: {0} Hours, {1} Minutes, {2} Seconds "
                            "#####################".format(hours, minutes, seconds))




