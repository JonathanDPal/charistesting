from parameter_test_infrastructure import *
import warnings
from time import time
from numpy import floor

start = time()
warnings.simplefilter('ignore', category=RuntimeWarning)

# Describe Fake Planets To Be Injected
fake_fluxes = [5e-3, 1e-3, 5e-4, 1e-4, 5e-5]
fake_seps = [20, 30, 40, 50, 60]

# KLIP Parameters To Be Sampled
annuli = [7,8,9]
subsections = [3,4,5]
movement = [1,2]
numbasis = [10,20,30]

# Filelist(s) & Associated Mask(s)
hd1160_fileset = 'HD1160_cubes/*.fits'
hd1160_mask = [144, 80]

# Create TestDataset For Each Set of Observations, e.g. hd = TestDataset('hd1160/*.fits', ...)
hd1160 = TestDataset(fileset=hd1160_fileset, object_name='HD1160', mask_xy=hd1160_mask,
                 fake_fluxes=fake_fluxes,fake_seps=fake_seps, annuli=annuli,
                 subsections=subsections, movement=movement, numbasis=numbasis)

# Have TestDataset Go To Town
hd1160.inject_fakes()
hd1160.run_KLIP()
hd1160.contrast_and_detection()

# Print Out Time Taken
end = time()
time_elapsed = end - start
hours = int(floor(time_elapsed / 3600))
remaining_time = time_elapsed - (hours * 3600)
minutes = int(floor(remaining_time / 60))
seconds = round(remaining_time - minutes * 60)

print("##################### TIME ELAPSED: {0} Hours, {1} Minutes, "
      "{2} Seconds #####################".format(hours, minutes, seconds))