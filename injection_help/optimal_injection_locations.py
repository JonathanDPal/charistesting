import numpy as np
from scipy.optimize import minimize, Bounds
from check_injection_locations import define_subsection_bounds, FWHMIOWA_calculator
from astropy.io import fits
import pandas as pd
from pyklip.klip import define_annuli_bounds

num_annuli = []  # a list of all values for annuli intended to be used
num_subsections = []  # a list of all values for subsections intended to be used
science_target_locs = []  # in form [x, y] or [[x1,y1], [x2,y2], ..., [xn, yn]]
file_with_sat_spots = ''  # a extracted datacube (NOT a KLIP output) with satellite spots written in to header
bounds_initialvals = []  # in form [(lb, ub, sep0), (lb, ub, pa0), ..., (lb, ub, sepn), (lb, ub, pan)], where n is
# equal to the number of fake planets you wish to inject (so this list should have 2n tuples of 3 elements)
output_file = ''

with fits.open(file_with_sat_spots) as hdulist:
    sat_spot_locs = [hdulist[1].header[hdr] for hdr in hdulist[1].header.keys() if 'SATS' in hdr]
    _, iwa, owa = FWHMIOWA_calculator(speccubefile=hdulist)
    SSxs, SSys = list(), list()
    for i, loc in enumerate(sat_spot_locs):
        parts = loc.split(' ')
        x, y = [float(part) for part in parts if part != '']
        SSxs.append(x)
        SSys.append(y)
    sat_spot_locs = [[x, y] for x, y in zip(SSxs, SSys)]

LB = np.array([elm[0] for elm in bounds_initialvals])
UB = np.array([elm[1] for elm in bounds_initialvals])
x0 = np.array([elm[2] for elm in bounds_initialvals])
bounds = Bounds(LB, UB)


def negdistance(fakelocs, annuli, subsections, st_locs, ss_locs, IWA, OWA):
    seps = [elm for idx, elm in enumerate(fakelocs) if idx % 2 == 0]
    pas = [elm for idx, elm in enumerate(fakelocs) if idx % 2 == 1]
    ann_boundaries = {ann: define_annuli_bounds(ann, IWA, OWA) for ann in annuli}
    sbs_boundaries = {sbs: define_subsection_bounds(sbs) for sbs in subsections}
    inj_locations = lsts_pasep_to_xy(pas, seps)

    xs = [loc[0] for loc in inj_locations]
    ys = [loc[1] for loc in inj_locations]

    if not isinstance(st_locs[0], (list, tuple, np.ndarray)):
        st_locs = [st_locs]

    for loc in st_locs:
        xs.append(loc[0])
        ys.append(loc[1])

    for loc in ss_locs:
        xs.append(loc[0])
        ys.append(loc[1])

    df = pd.DataFrame({'xpos': xs, 'ypos': ys})
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
            if all(row0) == all(row1):
                continue
            else:
                x0 = row0['xpos']
                y0 = row0['ypos']
                x1 = row1['xpos']
                y1 = row1['ypos']
                distances.append(np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2))
        min_distance_column.append(np.min(distances))

        anndistances = {ann: np.min(np.abs(np.array(ann_boundaries[ann]) - row0['sep'])) for ann in num_annuli}
        for ann in anndistances.keys():
            min_distance_from_ann[ann].append(anndistances[ann])
        sbsdistances = {sbs: row0['sep'] * np.min(np.abs(np.sin(np.array(sbs_boundaries[sbs]) - row0['pa']))) for sbs in
                        num_subsections}
        for sbs in sbsdistances.keys():
            min_distance_from_sbs[sbs].append(sbsdistances[sbs])

    everythingtomin = [np.min(min_distance_column)] + [np.min(lst) for lst in min_distance_from_ann.values()] + \
                      [np.min(lst) for lst in min_distance_from_sbs.values()]

    return -1 * np.min(everythingtomin)


result = minimize(fun=negdistance, args=(num_annuli, num_subsections, science_target_locs, sat_spot_locs, iwa, owa),
                  x0=x0, bounds=bounds)
if result.success:
    solution = result.x
    with open(output_file, 'w') as f:
        f.write('Sep, PA\n')
        seps = [elm for idx, elm in enumerate(solution) if idx % 2 == 0]
        pas = [elm for idx, elm in enumerate(solution) if idx % 2 == 1]
        for sep, pa in zip(seps, pas):
            f.write(f'{sep}, {pa}\n')
else:
    print(result.message)
