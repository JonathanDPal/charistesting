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
fileset0 =  # string to be passed into glob
mask0 =  # location(s) of science targets in [x-pos, y-pos]
object_name0 =  # string

# Setting Up Lists/Tuples For KLIP
annuli =  # list of integers
subsections =  # list of integers
movement =  # list of floats
spectrum =  # either None or 'methane'
numbasis =  # list of integers
corr_smooth =  # list of floats
highpass =  # list of floats or booleans

# Setting Mode For KLIP
mode =  # One of the following: 'ADI', 'SDI', 'ADI+SDI'

# Maximum Number of Threads to Use (set to None to use maximum computer capacity)
max_numthreads =  # integer or None

# Setting Up For Fake Planets
fake_fluxes =  # List(s) of Floats
fake_seps =  # List(s) of Integers
fake_PAs =  # List(s) of Integers

# Specifying Which Things to Do/Not Do (set to either True or False) #
# Most of the time, the four values below should be set to True
put_in_fakes =
run_KLIP_on_dataset_with_fakes =   # if no fakes are injected, this will just be a dataset without fakes
get_contrast =   # won't be calibrated if no fake planets are injected
get_planet_detections_from_dataset_with_fakes =
# Most of the time, these three values below should be set to False
run_KLIP_on_dataset_without_fakes =
get_planet_detections_from_dataset_without_fakes =
overwrite =  # whether or not to replace existing files if they exist

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
