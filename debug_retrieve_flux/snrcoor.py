from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import time
import sys

T = time.time()


def fit_func(X, a):
    return a * X


def valuefinder(filename):
    """
    Looks at a filename and discerns the KLIP parameters used to produce it. Can either find a specific KLIP
    parameter and return it in the form of a string, or it can find all KLIP parameters and return them in original
    form (int/float/bool).
    """
    starting_index = None
    final_index = None
    for h in range(len(filename)):
        if starting_index is None and filename[h] == '[':
            starting_index = h + 1
        elif final_index is None and filename[h] == ']':
            final_index = h
    return filename[starting_index: final_index].split(',')


def detections_filename(retrievals_filename):
    a, s, m, t, n, c, h = valuefinder(retrievals_filename)
    return f'HD1160/detections/{a}Annuli_{s}Subsections_{m}Movement_{t}Spectrum_{c}Smooth_{h}Highpass_-KL{n}_SNR-2.csv'


def distance(r, location):
    return np.sqrt((r['x'] - location[0]) ** 2 + (r['y'] - location[1]) ** 2)


files = glob('retrievals/*.xlsx')
# IGNORE THIS FOR NOW
if len(sys.argv) > 1:
    numfiles = len(files)
    batchsize = np.floor(float(numfiles) / 65.)
    start = sys.argv[1] * batchsize
    try:
        files = files[start: start + batchsize]
    except IndexError:
        files = files[start:]

vals = []
for i, file in enumerate(files):
    if i % 500 == 0:
        print(i)
    flux_x_y_snr = []
    df = pd.read_excel(file)
    for j, row in df.iterrows():
        if j > 5:
            break
        retfluxes = [row[col] for col in df.columns if 'Retrieved Flux' in col]
        xposes = [row[col] for col in df.columns if 'Retrieved X' in col]
        yposes = [row[col] for col in df.columns if 'Retrieved Y' in col]
        for rf, x, y, in zip(retfluxes, xposes, yposes):
            flux_x_y_snr.append([rf, x, y])

    snrfile = detections_filename(file)
    df1 = pd.read_csv(snrfile)
    fakes = df1[df1['Injected'] == 'True']
    xs = [elm[1] for elm in flux_x_y_snr]
    ys = [elm[2] for elm in flux_x_y_snr]
    for _, row in fakes.iterrows():
        found = False
        for k, loc in enumerate(zip(xs, ys)):
            if distance(row, loc) < 6.02:
                flux_x_y_snr[k].append(row['SNR Value'])
                found = True
                break
        if not found:
            flux_x_y_snr[i].append(0)
    for elm in flux_x_y_snr:
        vals.append(elm)

ret_fluxes = [elm[0] for elm in vals]
snr_vals = [elm[3] for elm in vals]

plt.scatter(ret_fluxes, snr_vals)
plt.savefig('snrcoor.png')
