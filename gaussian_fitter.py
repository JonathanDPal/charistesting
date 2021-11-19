import numpy as np
import sys
from astropy.io import fits

searchradius = 10  # optimize this based on how big you think a PSF should be

locs_file = sys.argv[1]
with open(locs_file) as f:
    lines = [line[:-1].replace(' ', '') if line[-1] == '\n' else line.replace(' ', '') for line in f]

img_file = sys.argv[2]
with fits.open(img_file) as f:
    data = f[0].data   # 1024 x 1024
data[np.where(np.isnan(data))] = 0

locs = [[int(z) for z in line.split(',')] for line in lines]

for i, loc in enumerate(locs):
    x0, y0 = loc
    searchbox = data[x0 - searchradius: x0 + searchradius + 1, y0 - searchradius: y0 + searchradius + 1]
    # center of PSF should be about at max value
    estimated_center = np.unravel_index(np.argmax(searchbox, axis=None), searchbox.shape)

    # building out arrays so that we get a full coordinate system when they are zipped together
    uppery, upperx = searchbox.shape
    yvals, xvals = np.arange(0, uppery, 1.0) - searchradius, np.arange(0, upperx, 1.0) - searchradius
    numyvals, numxvals = len(yvals), len(xvals)
    yvals = np.array(list(yvals) * numxvals)
    xvals = np.array([[xval] * numyvals for xval in xvals]).flatten()

    vals = [searchbox[int(y + searchradius), int(x + searchradius)] for y, x in zip(yvals, xvals)]
    vals_w_distances = pd.DataFrame({'vals': vals, 'distance squared': [y ** 2 + x ** 2 for y, x in zip(yvals, xvals)],
                                     'y': yvals, 'x': xvals})
    avg_vals_at_distance = {d2: np.mean(vals_w_distances[np.abs(vals_w_distances['distance squared']) == d2]) for d2
                            in np.unique(vals_w_distances['distance squared'])}
    peak = searchbox[estimated_center]
    below_half = [d2 if avg_vals_at_distance[d2] <= 0.5 * peak for d2 in avg_vals_at_distance.keys()]
    above_half = [d2 if avg_vals_at_distance[d2] >= 0.5 * peak for d2 in avg_vals_at_distance.keys()]
    if len(set(below_half).intersection(above_half)) > 0:
        fd2 = np.mean(list(set(below_half).intersection(above_half)))

