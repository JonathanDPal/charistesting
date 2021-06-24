from parameter_test_infrastructure import *
import warnings
from time import time
from numpy import floor
from astropy.io import fits
from glob import glob

# Keeping Track of Time Taken
start = time()

# KLIP often gets a bunch of RuntimeWarnings that we don't need to worry about
warnings.simplefilter('ignore', category=RuntimeWarning)

# Filelist(s) & Associated Mask(s) & Associated Name(s)
fileset0 = 'Fresh_HD1160_cubes/*.fits'
mask0 = None
object_name0 = 'FRESH'


# Describe Fake Planets To Be Injected
fake_fluxes = [1e-4, 1e-5, 1e-6]
fake_seps = [20, 40, 60]
fake_PAs=[0, 90, 180, 270]

with fits.open(glob(fileset0)[0]) as hdulist:
    filtname = hdulist[0].header['FILTNAME']
fake_fwhm = FWHMIOWA_calculator(speccubefile=None, filtname=filtname)[0]


# KLIP Parameters To Be Sampled
# annuli = [4, 6, 8, 10, 12]
annuli = [6]
# subsections = [2, 4, 6]
subsections = [2]
# movement = [0, 1, 2]
movement = [0]
# numbasis = [10, 20, 30, 40, 50, 60]
numbasis = [20]
# corr_smooth = [0, 1, 2]
corr_smooth = [1]
highpass = [True]
# highpass = [False, 5.0, True, 15.0] # True yields default (10.0)
spectrum = [None]
mode = 'ADI+SDI'

# Create TestDataset For Each Set of Observations
td0 = TestDataset(fileset=fileset0, object_name=object_name0, mask_xy=mask0, fake_fluxes=fake_fluxes,
                  fake_seps=fake_seps, annuli=annuli, subsections=subsections, movement=movement, numbasis=numbasis,
                  corr_smooth=corr_smooth, highpass=highpass, spectrum=spectrum, mode=mode, fake_PAs=fake_PAs,
                  fake_fwhm=fake_fwhm)

# Have TestDataset Run Each Part
td0.inject_fakes()
td0.run_KLIP(run_on_nofakes=False)
td0.contrast_and_detection(calibrate=[True])

# Print Out Time Taken
end = time()
time_elapsed = end - start
hours = int(floor(time_elapsed / 3600))
remaining_time = time_elapsed - (hours * 3600)
minutes = int(floor(remaining_time / 60))
seconds = round(remaining_time - minutes * 60)

print("##################### TIME ELAPSED: {0} Hours, {1} Minutes, {2} Seconds #####################".format(hours,
                                                                                                             minutes,
                                                                                                             seconds))
