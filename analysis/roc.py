import sys
from plotting_funcs import roc_generator, paramvaluesfinder
import numpy as np
import os

try:
    assert len(sys.argv) == 2  # the first argument in sys.argv is the name of the script ("roc.py")
except AssertionError:
    raise ValueError('Incorrect number of arguments. Should have one argument: either "highpass", "smooth", '
                     '"numbasis", or "all"')

param = sys.argv[1]  # parameter to look at (Highpass, Smooth, or Numbasis)

snr_vals = list(np.arange(start=3, stop=10.5, step=0.25))

# Using Log Files to Get Values
highpass_vals = paramvaluesfinder('Highpass')
smooth_vals = paramvaluesfinder('Corr_Smooth')
numbasis_vals = paramvaluesfinder('Numbasis')
ni = paramvaluesfinder('ni')  # number of injected planets

# Putting in to Tuples to Pass In to roc_generator
highpass = ('Highpass', highpass_vals)
smooth = ('Smooth', smooth_vals)
numbasis = ('KL', numbasis_vals)

# Output Filepaths
of1 = 'ROC/highpass_roc.png'
of11 = 'ROC/highpass_table.csv'
of2 = 'ROC/smooth_roc.png'
of22 = 'ROC/smooth_table.csv'
of3 = 'ROC/numbasis_roc.png'
of33 = 'ROC/numbasis_table.csv'

if not os.path.exists('ROC'):
    try:
        os.mkdir('ROC')
    except FileExistsError:  # if using analysis.sh, then very possible that two (simultaneous) scripts try to do
        # this at almost the exact same time
        pass

if str.lower(param) == 'highpass' or str.lower(param) == 'all':
    roc_generator(snr_vals, highpass, ni, of1)
    roc_generator(snr_vals, highpass, ni, of11, generate='table')
elif str.lower(param) == 'smooth' or str.lower(param) == 'all':
    roc_generator(snr_vals, smooth, ni, of2)
    roc_generator(snr_vals, smooth, ni, of22, generate='table')
elif str.lower(param) == 'numbasis' or str.lower(param) == 'all':
    roc_generator(snr_vals, numbasis, ni, of3)
    roc_generator(snr_vals, numbasis, ni, of33, generate='table')
