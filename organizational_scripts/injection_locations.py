import numpy as np
from astropy.io import fits
import pandas as pd


pas = []
seps = []
science_target_locs = []  # in form [x, y] or [[x1,y1], [x2,y2], ...]
file_with_sat_spots = ''  # a extracted datacube (NOT a KLIP output) with satellite spots written in to header


def pasep_to_xy(PAs, seps):
    """
    Takes lists of position angles and seperations and yields a numpy array with x-y coordinates for each combo, using
    the convention used in the table of values outputted by the planet detection software. The origin used in their
    convention is the center of the star's PSF.
    """
    PAs = [float(pa) for pa in PAs]  # if not a float, then pa / 180 will yield zero in many cases
    radians = np.array(PAs) / 180 * np.pi
    locs = []
    for sep in seps:
        for rad in radians:
            x = -np.sin(rad) * sep
            y = np.cos(rad) * sep
            loc = [x, y]
            locs.append(loc)
    return locs


def xy_to_pasep(xs, ys):
    seps = [np.sqrt(x ** 2 + y ** 2) for x, y in zip(xs, ys)]
    pas = [np.mod(np.rad2deg(np.arctan2(-x, y)), 360) for x, y in zip(xs, ys)]
    return [[sep, pa] for sep, pa in zip(seps, pas)]


def FWHMIOWA_calculator(speccubefile, filtname=None):
    """
    Finds FWHM, IWA, and OWA for a opened CHARIS data cube.
    """
    wavelengths = {'j': 1200e-9, 'h': 1550e-9, 'k': 2346e-9, 'broadband': 1550e-9}
    if filtname is None:
        wavelength = wavelengths[str.lower(speccubefile[1].header['FILTNAME'])]
    else:
        wavelength = wavelengths[str.lower(filtname)]
    D = 8
    lenslet_scale = 0.0162
    field_radius = 1.035
    FWHM = 2 * 1.22 * wavelength / D * 206265 / lenslet_scale
    IWA = 5
    OWA = (field_radius / lenslet_scale) - FWHM

    return FWHM, IWA, OWA


inj_locations = pasep_to_xy(pas, seps)

xs = [loc[0] for loc in inj_locations]
ys = [loc[1] for loc in inj_locations]

if not isinstance(science_target_locs[0], (list, tuple, np.ndarray)):
    science_target_locs = [science_target_locs]

xs.append(loc[0] for loc in science_target_locs)
ys.append(loc[1] for loc in science_target_locs)

with fits.open(file_with_sat_spots) as hdulist:
    sat_spot_locs = [hdulist[1].header[hdr] for hdr in hdulist[1].header.keys() if 'SATS' in hdr]
    dataset_fwhm = FWHMIOWA_calculator(hdulist)[0]

for i, loc in enumerate(sat_spot_locs):
    parts = loc.split(' ')
    x, y = [float(part) for part in parts if part != '']
    xs.append(x)
    ys.append(y)

types = []
for _ in range(len(inj_locations)):
    types.append('fake planet')
for _ in range(len(science_target_locs)):
    types.append('science target')
for _ in range(len(sat_spot_locs)):
    types.append('satellite spot')

df = pd.DataFrame((xs, ys, types), columns=['xpos', 'ypos', 'type'])

min_distance_sq_column = []
for _, row0 in df.iterrows():
    distances2 = []
    for _, row1 in df.iterrows():
        if row0 == row1:
            continue
        else:
            x0 = row0['xpos']
            y0 = row0['ypos']
            x1 = row1['xpos']
            y1 = row1['ypos']
            distances2.append((x1 - x0) ** 2 + (y1 - y0) ** 2)
    min_distance_sq_column.append(np.min(distances2))

df['closest distance squared'] = min_distance_sq_column

if np.min(min_distance_sq_column) < 2 * dataset_fwhm:
    print("Some things are too close together. See information below")
    print(df)
else:
    print("Everything is at least 2 FWHMs apart!")
