import sys
from glob import glob
import os

direc = sys.argv[1]
klipfiles = glob(f'{direc}/klipped_cubes_W_fakes/*.fits')
contrastfiles = glob(f'{direc}/calibrated_contrast/*.csv')
detectionsfiles = glob(f'{direc}/detections/*.csv')

files = klipfiles + contrastfiles + detectionsfiles
for file in files:
    os.system(f'rsync -caqz {file} jpal@planetfinder.esc.nd.edu:/home/jpal/data0/jpal/parameter_sampling/'
              f'{direc}/from_crc')
    os.remove(file)