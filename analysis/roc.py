import sys
from plotting_funcs import roc_generator, paramvaluesfinder
import numpy as np
import os

try:
    assert len(sys.argv) == 4
except AssertionError:
    raise ValueError("Incorrect number of arguments. First argument should be either 'highpass', 'smooth', "
                     "or 'numbasis', and second argument should be the number of planets which were injected.")

name = sys.argv[1]  # name of directory with all the data in it
param = sys.argv[2]  # parameter to look at (Highpass,
ni = int(sys.argv[3])  # number of injected planets

snr_vals = list(np.arange(start=3, stop=10.5, step=0.25))

# Using Log Files to Get Values
highpass_vals = paramvaluesfinder(name, 'Highpass')
smooth_vals = paramvaluesfinder(name, 'Corr_Smooth')
numbasis_vals = paramvaluesfinder(name, 'Numbasis')

# Putting in to Tuples to Pass In to roc_generator
highpass = ("Highpass", highpass_vals)
smooth = ("Smooth", smooth_vals)
numbasis = ("KL", numbasis_vals)

# Output Filepaths
of1 = 'ROC/highpass_roc.png'
of2 = 'ROC/smooth_roc.png'
of3 = 'ROC/numbasis_roc.png'

if not os.path.exists('ROC'):
    os.mkdir('ROC')

if param == 'highpass':
    roc_generator(snr_vals, highpass, ni, name, of1)
elif param == 'smooth':
    roc_generator(snr_vals, smooth, ni, name, of2)
elif param == 'numbasis':
    roc_generator(snr_vals, numbasis, ni, name, of3)
