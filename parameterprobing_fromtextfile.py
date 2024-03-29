import numpy as np
from parameter_test_infrastructure import *
from time import time
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

# Text File to Get Params From
text_file =

# Setting Mode For KLIP
mode =  # One of the following: 'ADI', 'SDI', 'ADI+SDI'

# Maximum Number of Threads to Use (set to None to use maximum computer capacity)
max_numthreads =  # integer or None

# If True, we'll use the version of pyKLIP which uses less memory
memorylite =  # True or False

# Setting Up For Fake Planets
fake_fluxes =  # List(s) of Floats
fake_seps =  # List(s) of Integers
fake_PAs =  # List(s) of Integers
compiler =  # either 'zip' or 'iterate'
numsepgroups =  # integer if compiler == 'zip'; None if compiler == 'iterate'
tweak_injections =  # boolean

# Specifying Which Things to Do/Not Do (set to either True or False) #
# Most of the time, the five values below should be set to True
put_in_fakes =
run_KLIP_on_dataset_with_fakes =   # if no fakes are injected, this will just be a dataset without fakes
get_contrast =   # won't be calibrated if no fake planets are injected
get_planet_detections_from_dataset_with_fakes =
create_log_file =  # whether or not to create log file; note that if one exists, then it will be destroyed and replaced
# Most of the time, these four values below should be set to False
run_KLIP_on_dataset_without_fakes =
get_planet_detections_from_dataset_without_fakes =
airydisk_kernel = # for cross-correlation prior to creating SNR map for point source detections (if False, then Gauss)
overwrite =   # whether or not to replace existing files if they exist
verbose =  # if False, then less stuff printed and less stuff in log file

######################
# END OF USER INPUTS #
######################

########################
# CHECKING USER INPUTS #
########################

# Reading From Text File to Get Parameters
annuli, subsections, movement, spectrum, numbasis, corr_smooth, highpass = params_from_text_file(text_file)

# Converting into list format if single value given and sanitizing format
if isinstance(annuli, (int, float)):
    annuli = [int(annuli)]
if isinstance(subsections, (int, float)):
    subsections = [int(subsections)]
if isinstance(movement, (int, float)):
    movement = [float(movement)]
if isinstance(spectrum, (type(None), str)):
    spectrum = [spectrum]
if isinstance(numbasis, (int, float)):
    numbasis = [int(numbasis)]
if isinstance(corr_smooth, (int, float)):
    corr_smooth = [float(corr_smooth)]
if isinstance(highpass, (bool, int, float)):
    highpass = [float(highpass) if not isinstance(highpass, bool) else highpass]

for param in [[annuli, 'annuli'], [subsections, 'subsections'], [movement, 'movement'], [spectrum, 'spectrum'],
              [corr_smooth, 'corr_smooth'], [highpass, 'highpass']]:
    if not isinstance(param[0], (list, tuple, np.ndarray)):
        raise TypeError(f"{param[1]} could not be coerced into a list safely. Check input. See "
                        "https://docs.google.com/document/d/1yX0l96IZs1IxxKCRmriVSAQM3KFGF9U1-FnpJXhcLXo/edit?usp"
                        "=sharing for help")

# Common Mistake is Inputting Mode as a List/Tuple Like Other Params
if not isinstance(mode, str):
    raise TypeError("Mode needs to be a string. Check input.")

if overwrite not in [True, False]:
    raise TypeError("Overwrite needs to be either True or False. Check input.")

# Compressing info about fake planet injections into desired format
if compiler == 'zip':
    fakes = [np.array((flux, sep, pa)) for flux, sep, pa in zip(fake_fluxes, fake_seps, fake_PAs)]
elif compiler == 'iterate':
    fakes = list()
    if type(fake_fluxes[0]) not in [tuple, list]:
        fake_fluxes = [fake_fluxes]
    if type(fake_PAs[0]) not in [tuple, list]:
        fake_PAs = [fake_PAs]
    for fluxgroup, pagroup in zip(fake_fluxes, fake_PAs):
        for fake_flux, sep in zip(fluxgroup, fake_seps):
            for pa in pagroup:
                fakes.append(np.array((fake_flux, sep, pa)))
    numsepgroups = len(fake_seps)
elif compiler is None:
    fakes = None
else:
    raise ValueError('Check value for "compiler" variable.')

if fakes is not None:
    assert len(fakes) % numsepgroups == 0, 'Check your fake planet specifications. The number of seperation groups ' \
                                           'does not divide the number of fake planets.'

# SYNTHESIZING USER INPUTS INTO A COUPLE ADDITIONAL THINGS #
detect_planets = get_planet_detections_from_dataset_with_fakes or get_planet_detections_from_dataset_without_fakes
get_contrast_and_detections = get_contrast or detect_planets
if run_KLIP_on_dataset_without_fakes or run_KLIP_on_dataset_with_fakes:
    build_charis_data = 'true'
else:
    build_charis_data = 'temporary'
if airydisk_kernel:
    kernel_type = 'airy'
else:
    kernel_type = 'gaussian'

# Selecting Batch if Needed (or providing link to instructions document)
if len(sys.argv) != 1:
    if "-h" in sys.argv or "--help" in sys.argv or "h" in sys.argv or "help" in sys.argv:
        raise KeyboardInterrupt("See https://docs.google.com/document/d/1yX0l96IZs1IxxKCRmriVSAQM3KFGF9U1-FnpJXhcLXo"
                                 "/edit?usp=sharing for help.")
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
fake_fwhm0 = FWHMIOWA_calculator(glob(fileset0))[0]

# Create TestDataset For Target 0
td0 = TestDataset(fileset=fileset0, object_name=object_name0, mask_xy=mask0, fakes=fakes, numsepgroups=numsepgroups,
                  annuli=annuli, subsections=subsections, movement=movement, numbasis=numbasis, corr_smooth=corr_smooth,
                  highpass=highpass, spectrum=spectrum, mode=mode, fake_fwhm=fake_fwhm0, batched=batched,
                  overwrite=overwrite, memorylite=memorylite, build_all_combos=False,
                  build_charis_data=build_charis_data, verbose=verbose, generatelogfile=create_log_file,
                  tweak_injections=tweak_injections)

# Have TestDataset 0 Run Each Part
# if we want KLIP output of data without fakes, we need to run KLIP before injecting planets
if run_KLIP_on_dataset_without_fakes:
    td0.run_KLIP_on_data_without_fakes(numthreads=max_numthreads)
if run_KLIP_on_dataset_with_fakes:
    td0.inject_fakes()
    td0.run_KLIP_on_data_with_fakes(numthreads=max_numthreads)
if get_contrast_and_detections:
    td0.contrast_and_detection(run_contrast=get_contrast, run_planet_detection=detect_planets,
                               datasetwithfakes=put_in_fakes, kernel_type=kernel_type)

# Print Out Time Taken
end0 = time()
time_elapsed = end0 - start0
hours = int(np.floor(time_elapsed / 3600))
remaining_time = time_elapsed - (hours * 3600)
minutes = int(np.floor(remaining_time / 60))
seconds = round(remaining_time - minutes * 60)

if create_log_file and verbose:
    td0.write_to_log_and_print(f'##################### TIME ELAPSED: {hours} Hours, {minutes} Minutes, {seconds} Seconds '
                                '#####################')
elif verbose:
    print(f'##################### TIME ELAPSED: {hours} Hours, {minutes} Minutes, {seconds} Seconds '
          f'#####################')
if not create_log_file:
    os.remove(f'{object_name0}/temp.txt')