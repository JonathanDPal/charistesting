from parameter_test_infrastructure import *
from time import time
from numpy import floor
from astropy.io import fits
from glob import glob
import os
import sys
import warnings
from astropy.wcs.wcs import FITSFixedWarning

###############
# USER INPUTS #
###############

# see https://docs.google.com/document/d/1yX0l96IZs1IxxKCRmriVSAQM3KFGF9U1-FnpJXhcLXo/edit?usp=sharing for help

# General Set-Up
fileset0 = 'HD1160_cubes/*.fits'  # string to be passed into glob
mask0 = [144, 80]  # location(s) of science targets in [x-pos, y-pos]
object_name0 = 'HD1160'  # string

fileset1 = 'HR8799_cubes/*.fits'
mask1 = [[128, 152], [157, 137], [129, 70]]
object_name1 = 'HR8799'

# Setting Up Lists/Tuples For KLIP
annuli = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]  # list of integers
subsections = [2, 4, 6]  # list of integers
movement = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]  # list of floats
spectrum = [None, 'methane']  # list containing None and/or 'methane'
numbasis = [10, 20, 30, 40, 50, 60]  # list of integers
corr_smooth = [0.0, 1.0, 2.0, 3.0]  # list of floats
highpass = [False, 5.0, True]  # list of floats or booleans (note: True & 10.0 are equivalent)

# Setting Mode For KLIP
mode = 'ADI+SDI'  # One of the following: 'ADI', 'SDI', 'ADI+SDI'

# Maximum Number of Threads to Use (set to None to use maximum computer capacity)
max_numthreads = None  # integer or None

# Setting Up For Fake Planets
# in this example, there is one tier of planets, 4 injected at each seperation (12 total fake planets)
fake_fluxes = [1e-4, 1e-5, 1e-6]  # List(s) of Floats
fake_seps = [20, 40, 60]  # List(s) of Integers
fake_PAs = [0, 90, 180, 270]  # List(s) of Integers

# Specifying Which Things to Do/Not Do (set to either True or False) #
# Most of the time, the four values below should be set to True
put_in_fakes = True
run_KLIP_on_dataset_with_fakes = True  # if no fakes are injected, this will just be a dataset without fakes
get_contrast = True  # won't be calibrated if no fake planets are injected
get_planet_detections_from_dataset_with_fakes = True
# Most of the time, these three values below should be set to False
run_KLIP_on_dataset_without_fakes = False
get_planet_detections_from_dataset_without_fakes = False
overwrite = False  # whether or not to replace existing files if they exist

######################
# END OF USER INPUTS #
######################

########################
# CHECKING USER INPUTS #
########################

# Making Sure This Group of Parameters Are In The Form of a List
for param in [[annuli, 'annuli'], [subsections, 'subsections'], [movement, 'movement'], [spectrum, 'spectrum'],
              [corr_smooth, 'corr_smooth'], [highpass, 'highpass']]:
    if not isinstance(param[0], (list, tuple)):
        raise TypeError(f"{param[1]} needs to be a list. Check input. See "
                        "https://docs.google.com/document/d/1yX0l96IZs1IxxKCRmriVSAQM3KFGF9U1-FnpJXhcLXo/edit?usp"
                        "=sharing for help")

# Checking Mode -- Common Mistake is Inputting Mode as a List/Tuple Like Other Params
if not isinstance(mode, str):
    raise TypeError("Mode needs to be a string. Check input. See "
                    "https://docs.google.com/document/d/1yX0l96IZs1IxxKCRmriVSAQM3KFGF9U1-FnpJXhcLXo/edit?usp"
                    "=sharing for help")

# SYNTHESIZING USER INPUTS INTO A COUPLE ADDITIONAL BOOLEANS #
detect_planets = get_planet_detections_from_dataset_with_fakes or get_planet_detections_from_dataset_without_fakes
if put_in_fakes:
    datasetwithfakes = True
else:
    datasetwithfakes = False
get_contrast_and_detections = get_contrast or detect_planets

# Selecting Batch if Needed (or providing link to instructions document)
if len(sys.argv) != 1:
    if "-h" in sys.argv or "--help" in sys.argv:
        print("See https://docs.google.com/document/d/1yX0l96IZs1IxxKCRmriVSAQM3KFGF9U1-FnpJXhcLXo/edit?usp=sharing "
              "for help.")
        raise KeyboardInterrupt
    elif len(sys.argv) == 3:
        batchindex = int(sys.argv[1])
        batchsize = int(sys.argv[2])
        batched = (True, batchindex, batchsize)
    else:
        raise ValueError("Incorrect number of keyword arguments passed in. Please either provide two (batchindex & "
                         "batchsize) or zero.")
else:
    batched = False

#############################
# ACTUALLY STARTING TESTING #
#############################
start0 = time()

# KLIP creates a bunch of RuntimeWarnings that we don't want to spam our log file
warnings.simplefilter('ignore', category=RuntimeWarning)
warnings.simplefilter('ignore', category=FITSFixedWarning)

# FWHM Is Based on Central Wavelength of Filter Used During Obsveration.
with fits.open(glob(fileset0)[0]) as hdulist:
    fake_fwhm0 = FWHMIOWA_calculator(hdulist)[0]

# Create TestDataset For Target 0
td0 = TestDataset(fileset=fileset0, object_name=object_name0, mask_xy=mask0, fake_fluxes=fake_fluxes,
                  fake_seps=fake_seps, annuli=annuli, subsections=subsections, movement=movement, numbasis=numbasis,
                  corr_smooth=corr_smooth, highpass=highpass, spectrum=spectrum, mode=mode, fake_PAs=fake_PAs,
                  fake_fwhm=fake_fwhm0, batched=batched, overwrite=overwrite)

# Have TestDataset 0 Run Each Part
# if we want KLIP output of data without fakes, we need to run KLIP before injecting planets
if run_KLIP_on_dataset_without_fakes:
    td0.run_KLIP_on_data_without_fakes(numthreads=max_numthreads)
if put_in_fakes:
    td0.inject_fakes()
if run_KLIP_on_dataset_with_fakes:
    td0.run_KLIP_on_data_with_fakes(numthreads=max_numthreads)
if get_contrast_and_detections:
    td0.contrast_and_detection(run_planet_detection=detect_planets, datasetwithfakes=datasetwithfakes,
                               numthreads=max_numthreads)

# Print Out Time Taken
end0 = time()
time_elapsed = end0 - start0
hours = int(floor(time_elapsed / 3600))
remaining_time = time_elapsed - (hours * 3600)
minutes = int(floor(remaining_time / 60))
seconds = round(remaining_time - minutes * 60)

td0.write_to_log_and_print(f'##################### TIME ELAPSED: {hours} Hours, {minutes} Minutes, {seconds} Seconds '
                           '#####################')

# Taking Stuff Out of Memory
td0.dataset = None
for i in range(len(td0.trials)):
    td0.trials[i] = None
td0 = None

# STARTING SECOND OBSERATION SET #
start1 = time()

with fits.open(glob(fileset1)[0]) as hdulist:
    fake_fwhm1 = FWHMIOWA_calculator(hdulist)[0]

td1 = TestDataset(fileset=fileset1, object_name=object_name1, mask_xy=mask1, fake_fluxes=fake_fluxes,
                  fake_seps=fake_seps, annuli=annuli, subsections=subsections, movement=movement, numbasis=numbasis,
                  corr_smooth=corr_smooth, highpass=highpass, spectrum=spectrum, mode=mode, fake_PAs=fake_PAs,
                  fake_fwhm=fake_fwhm1, batched=batched, overwrite=overwrite)

if run_KLIP_on_dataset_without_fakes:
    td1.run_KLIP_on_data_without_fakes(numthreads=max_numthreads)
if put_in_fakes:
    td1.inject_fakes()
if run_KLIP_on_dataset_with_fakes:
    td1.run_KLIP_on_data_with_fakes(numthreads=max_numthreads)
if get_contrast_and_detections:
    td1.contrast_and_detection(run_planet_detection=detect_planets, datasetwithfakes=datasetwithfakes,
                               numthreads=max_numthreads)

end1 = time()
time_elapsed = end1 - start1
hours = int(floor(time_elapsed / 3600))
remaining_time = time_elapsed - (hours * 3600)
minutes = int(floor(remaining_time / 60))
seconds = round(remaining_time - minutes * 60)

td1.write_to_log_and_print(f'##################### TIME ELAPSED: {hours} Hours, {minutes} Minutes, {seconds} Seconds '
                           '#####################')
