from parameter_test_infrastructure import *
import warnings
from time import time
from numpy import floor

start = time()
warnings.simplefilter('ignore', category=RuntimeWarning)

# Describe Fake Planets To Be Injected
fake_fluxes = [5e-3, 5e-4, 5e-5]
fake_seps = [20, 40, 60]

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
highpass = [False]
spectrum = [None]
mode = 'ADI+SDI'

# Filelist(s) & Associated Mask(s)
fileset = 'HR8799_cubes/*.fits'
mask = None
object_name = 'HR8799_mod2'

# Create TestDataset For Each Set of Observations, e.g. hd = TestDataset('hd1160/*.fits', ...)
td1 = TestDataset(fileset=fileset, object_name=object_name, mask_xy=mask,
                 fake_fluxes=fake_fluxes,fake_seps=fake_seps, annuli=annuli,
                 subsections=subsections, movement=movement, numbasis=numbasis,
                     corr_smooth=corr_smooth, highpass=highpass, spectrum=spectrum, mode=mode,
                  fake_PAs=[20,110,200,290])

# Have TestDataset Go To Town
td1.inject_fakes()
td1.run_KLIP()
td1.contrast_and_detection()

# Print Out Time Taken
end = time()
time_elapsed = end - start
hours = int(floor(time_elapsed / 3600))
remaining_time = time_elapsed - (hours * 3600)
minutes = int(floor(remaining_time / 60))
seconds = round(remaining_time - minutes * 60)

print("##################### TIME ELAPSED: {0} Hours, {1} Minutes, "
      "{2} Seconds #####################".format(hours, minutes, seconds))
