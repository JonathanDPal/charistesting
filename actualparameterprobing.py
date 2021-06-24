from parameter_test_infrastructure import *
import warnings
from time import time
from numpy import floor

# Keeping Track of Time Taken
start = time()

# KLIP often gets a bunch of RuntimeWarnings that we don't need to worry about
warnings.simplefilter('ignore', category=RuntimeWarning)

# Describe Fake Planets To Be Injected
fake_fluxes = [1e-4]
fake_seps = [40]
fake_fwhm = 6.01924560185 # broadband
fake_PAs=[0]

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
# highpass = [False, 5.0, True, 15.0]
spectrum = [None]
mode = 'ADI+SDI'

# Filelist(s) & Associated Mask(s) & Associated Name(s)
fileset0 = 'HR8799_cubes/*.fits'
mask0 = None
object_name0 = 'HR8799_1pt'

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

print(f"##################### TIME ELAPSED: {hours} Hours, {minutes} Minutes, {seconds} Seconds #####################")
