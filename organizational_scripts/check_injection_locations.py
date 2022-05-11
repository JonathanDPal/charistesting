import numpy as np
from astropy.io import fits
import pandas as pd
from pyklip.klip import define_annuli_bounds


num_annuli = []  # a list of all values for annuli intended to be used
num_subsections = []  # a list of all values for subsections intended to be used
# the script is going to assume that every combination of position angle and seperation is intended to be used
pas = []  # position angles where you intend to put in fake planets
seps = []  # seperations where you intend to put in fake planets
science_target_locs = []  # in form [x, y] or [[x1,y1], [x2,y2], ..., [xn, yn]]
file_with_sat_spots = ''  # a extracted datacube (NOT a KLIP output) with satellite spots written in to header


def define_subsection_bounds(sbs):
    """
    Taken verbatim from pyklip.parallelized.klip_dataset()
    """
    dphi = 2 * np.pi / subsections
    phi_bounds = [[dphi * phi_i - np.pi, dphi * (phi_i + 1) - np.pi] for phi_i in range(sbs)]
    phi_bounds[-1][1] = np.pi
    return phi_bounds


ann_boundaries = {ann: define_annuli_bounds(ann, 5, OWA) for ann in num_annuli}
sbs_boundaries = {sbs: define_subsection_bounds(sbs) for sbs in num_subsections}


def single_pasep_to_xy(rad, sep):
    x = -np.sin(rad) * sep
    y = np.cos(rad) * sep
    return x, y


def lsts_pasep_to_xy(PAs, seps):
    """
    Takes lists of position angles and seperations and yields a numpy array with x-y coordinates for each combo, using
    the convention used in the table of values outputted by the planet detection software. The origin used in their
    convention is the center of the star's PSF.
    """
    PAs = [float(pa) for pa in PAs]  # if not a float, then pa / 180 will yield zero in many cases
    radians = np.array(PAs) / 180 * np.pi
    locs = list()
    for sep in seps:
        for rad in radians:
            x, y = single_pasep_to_xy(rad, sep)
            loc = [x, y]
            locs.append(loc)
    return locs


def single_xy_to_pasep(x, y):
    sep = np.sqrt(x ** 2 + y ** 2)
    pa = np.mod(np.rad2deg(np.arctan2(-x, y)), 360)
    return sep, pa


def lsts_xy_to_pasep(xs, ys):
    seps = [np.sqrt(x ** 2 + y ** 2) for x, y in zip(xs, ys)]
    pas = [np.mod(np.rad2deg(np.arctan2(-x, y)), 360) for x, y in zip(xs, ys)]
    return [[sep, pa] for sep, pa in zip(seps, pas)]


def FWHMIOWA_calculator(speccubefile, filtname=None):
    """
    This used to calculate FWHM based on the central wavelength, but we found that the equation for that calculation
    was faulty, so at this point, it just spits back 3.5 as the FWHM. Just kept in this format to minimize the
    number of edits needed and also to permit modifications in the future such that we can have different FWHMs for
    different central wavelengths.
    """
    wavelengths = {'j': 1200e-9, 'h': 1550e-9, 'k': 2346e-9, 'broadband': 1550e-9}
    if filtname is None:
        wavelength = wavelengths[str.lower(speccubefile[1].header['FILTNAME'])]
    else:
        wavelength = wavelengths[str.lower(filtname)]
    D = 8
    # In the future, we'll probably have FWHM look something like this:
    # fwhms = {'j': {fwhm1}, 'h': {fwhm2}, 'k':{fwhm3}, 'broadband': 3.5}
    # FWHM = fwhms[wavelength]
    lenslet_scale = 0.0162
    field_radius = 1.035
    FWHM = 3.5
    IWA = 5
    OWA = (field_radius / lenslet_scale) - FWHM

    return FWHM, IWA, OWA


inj_locations = lsts_pasep_to_xy(pas, seps)

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

types = list()
for _ in range(len(inj_locations)):
    types.append('fake planet')
for _ in range(len(science_target_locs)):
    types.append('science target')
for _ in range(len(sat_spot_locs)):
    types.append('satellite spot')

df = pd.DataFrame((xs, ys, types), columns=['xpos', 'ypos', 'type'])
dfsepspas = lsts_xy_to_pasep(xs, ys)
dfseps = [elm[0] for elm in dfsepspas]
dfpas = [elm[1] for elm in dfsepspas]
df['sep'] = dfseps
df['pa'] = dfpas

min_distance_column = list()
min_distance_from_sbs = {sbs: list() for sbs in num_subsections}
min_distance_from_ann = {ann: list() for ann in num_annuli}
for _, row0 in df.iterrows():
    distances = list()
    for _, row1 in df.iterrows():
        if row0 == row1:
            continue
        else:
            x0 = row0['xpos']
            y0 = row0['ypos']
            x1 = row1['xpos']
            y1 = row1['ypos']
            distances.append(np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2))
    min_distance_column.append(np.min(distances))

    anndistances = {ann: np.min(np.abs(np.array(ann_boundaries[ann])) - row0['sep']) for ann in num_annuli}
    for ann in anndistances.keys():
        min_distance_from_ann[ann].append(anndistances[ann])
    sbsdistances = {sbs: row0['sep'] * np.min(np.abs(np.sin(np.array(sbs_boundaries[sbs]) - row0['pa']))) for sbs in
                    num_subsections}
    for sbs in sbsdistances.keys():
        min_distance_from_sbs[sbs].append(sbsdistances[sbs])

df['min distance from objects'] = min_distance_column
for ann in min_distance_from_ann.keys():
    df[f'min distance from radial bounds for ann={ann}'] = min_distance_from_ann[ann]
for sbs in min_distance_from_sbs.keys():
    df[f'min distance from angle bounds for sbs={sbs}'] = min_distance_from_sbs[sbs]

everythingtomin = [np.min(min_distance_column)] + [np.min(lst) for lst in min_distance_from_ann.values()] + \
                  [np.min(lst) for lst in min_distance_from_sbs.values()]
if np.min(everythingtomin) < 2 * dataset_fwhm:
    print(f"Some things are too close together. See information below, specifically the value(s) in the right most "
          f"columns which are less than {2 * dataset_fwhm}.")
    print(df)
else:
    print('Looks good!')
