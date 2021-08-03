from pyklip.klip import define_annuli_bounds
import numpy as np
from glob import glob
from pyklip.instruments.CHARIS import CHARISData
from astropy.wcs.wcs import FITSFixedWarning
import warnings
import sys

warnings.simplefilter('ignore', category=RuntimeWarning)
warnings.simplefilter('ignore', category=FITSFixedWarning)


def define_subsection_bounds(sbs):
    """
    Taken verbatim from pyklip.parallelized.klip_dataset()
    """
    dphi = 2 * np.pi / subsections
    phi_bounds = [[dphi * phi_i - np.pi, dphi * (phi_i + 1) - np.pi] for phi_i in range(sbs)]
    phi_bounds[-1][1] = np.pi
    return tuple(phi_bounds)


dataset = CHARISData(glob('cubes/*.fits'))
imgs = dataset.input
centers = dataset.centers

IWA = 5
OWA = None
annuli, subsections = int(sys.argv[1]), int(sys.argv[1])

dims = imgs.shape
x, y = np.meshgrid(np.arange(dims[2] * 1.0), np.arange(dims[1] * 1.0))
nanpix = np.where(np.isnan(imgs[0]))
# if user didn't supply how to define NaNs
if OWA is None:
    full_image = True  # reduce the full image
    # define OWA as either the closest NaN pixel or edge of image if no NaNs exist
    if np.size(nanpix) == 0:
        OWA = np.sqrt(np.max((x - centers[0][0]) ** 2 + (y - centers[0][1]) ** 2))
    else:
        # grab the NaN from the 1st percentile (this way we drop outliers)
        OWA = np.sqrt(np.percentile((x[nanpix] - centers[0][0]) ** 2 + (y[nanpix] - centers[0][1]) ** 2, 1))
else:
    full_image = False  # don't reduce the full image, only up the the IWA

rb = define_annuli_bounds(annuli, IWA=IWA, OWA=OWA)
if full_image:
    rb[annuli - 1] = np.nan

rad_bounds = []
for r in rb:
    if isinstance(r, (tuple, list)):
        for elm in r:
            rad_bounds.append(round(elm, 2))
    else:
        rad_bounds.append(round(r, 2))

sb = define_subsection_bounds(subsections)
ang_bounds = []
for s in sb:
    if isinstance(s, (tuple, list)):
        for elm in s:
            ang_bounds.append(round(elm, 2))
    else:
        ang_bounds.append(round(s, 2))
ang_bounds = list(set(ang_bounds))

with open('as.txt', 'a') as file:
    file.write(f'\nAnnuli = {annuli} -- {rad_bounds}')

with open('as.txt', 'a') as file:
    file.write(f'\nSubsections = {subsections} -- {ang_bounds}')
