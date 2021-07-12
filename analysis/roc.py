from plotting_funcs import roc_generator
import numpy as np

snr_vals = np.arange(start=3, stop=10.5, step=0.25)

highpass = ("Highpass", [False, 5.0, True, 15.0])
smooth = ("Smooth", [0, 1, 2])
numbasis = ("KL", [10,20,30,40,50,60])

name = "HD1160" 

of1 = 'ROC/highpass_roc.png'
of2 = 'ROC/smooth_roc.png'
of3 = 'ROC/numbasis_roc.png'

num_injections = 18

#roc_generator(snr_vals, highpass, num_injections, name, of1)
roc_generator(snr_vals, smooth, num_injections, name, of2)
#roc_generator(snr_vals, numbasis, num_injections, name, of3)

