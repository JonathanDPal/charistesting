from astropy.io import fits
import sys
from pyklip.instruments.CHARIS import CHARISData
from parameter_test_infrastructure import make_dn_per_contrast
from glob import glob

fileset0 = 'HIP86032_cubes/*.fits'
object_name0 = 'HIP86032'

dataset = make_dn_per_contrast(CHARISData(glob(fileset0)))
rot_angs = dataset.PAs
flipx = dataset.flipx
dn_per_contrast = dataset.dn_per_contrast
wln_um = dataset.wvs
length = dataset.input.shape[1]

with open(f'{object_name0}_dataset_info.txt', 'w') as f:
    f.write(f'Rotation Angles:\n')
    for r_ang in rot_angs:
        f.write(f'{r_ang}\n')
    f.write(f'Flip_x: {flipx}\n')
    f.write(f'DN_per_contrast:\n')
    for ratio in dn_per_contrast:
        f.write(f'{ratio}\n')
    f.write(f'Wavelengths:\n')
    for wln in wln_um:
        f.write(f'{wln}\n')
    f.write(f'Length: {length}')
