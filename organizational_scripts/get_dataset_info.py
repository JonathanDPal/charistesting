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
    f.write(f'Rotation Angles: {rot_angs}, {type(rot_angs)}\n')
    f.write(f'Flip_x: {flipx}, {type(flipx)}\n')
    f.write(f'DN_per_contrast: {dn_per_contrast}, {type(dn_per_contrast)}\n')
    f.write(f'Wavelengths: {wln_um}, {type(wln_um)}\n')
    f.write(f'Length: {length}, {type(length)}')
    for thing in [rot_angs, flipx, dn_per_contrast, wln_um, length]:
        if type(thing) == list:
            f.write(f'{type(thing[0])}\n')
