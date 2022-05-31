import sys
from glob import glob
import os
import itertools
import numpy as np

direc = sys.argv[1]
klipfiles = glob(f'{direc}/klipped_cubes_W_fakes/*.fits')
klipnofakesfiles = glob(f'{direc}/klipped_cubes_N_fakes/*.fits')
contrastfiles = glob(f'{direc}/calibrated_contrast/*.csv')
uncalcontrastfiles = glob(f'{direc}/uncalibrated_contrast/*.csv')
detectionsfiles = glob(f'{direc}/detections/*.csv')

initfiles = klipfiles + klipnofakesfiles + contrastfiles + uncalcontrastfiles + detectionsfiles
perm = list(itertools.permutations(np.arange(len(initfiles))))[np.random.randint(len(initfiles))]
newfiles = [initfiles[perm[idx]] for idx in range(len(initfiles))]
for file in newfiles:
    subdirec = file.split('/')[-2]
    if subdirec in ['klipped_cubes_W_fakes', 'klipped_cubes_N_fakes']:
        try:
            assert os.path.exists(file)
            os.system(f'rsync -caqz {file} jpal@planetfinder.esc.nd.edu:/home/jpal/data1/jpal/{direc}/{subdirec}')
            os.remove(file)
        except (AssertionError, FileNotFoundError) as e:
            pass
    else:
        try:
            assert os.path.exists(file)
            os.system(f'rsync -caqz {file} jpal@planetfinder.esc.nd.edu:/home/jpal/data0/jpal/{direc}/{subdirec}')
            os.remove(file)
        except (AssertionError, FileNotFoundError) as e:
            pass
