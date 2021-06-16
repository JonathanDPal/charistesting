import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from glob import glob

snr_values = np.arange(start=3, stop=10, step=0.5)

def roc_generator(snr_values, vals_and_finders, filepath_to_save):
    """
    Makes an ROC curve.
    ---
    Args:
        snr_values (list): SNR Values to plot
        vals_and_finders (dict): Values to check outare keys and file finders are the values.
        filepath_to_save (str): Where to save graph.
    """
    collection = dict()
    for val in vals_and_finders.keys():
        collection[val] = dict()
        for snr_value in snr_values:
            collection[val][snr_value] = list()
        val_fileset = glob(vals_and_finders[val])
        for file in val_fileset:
            detections = pd.read_csv(file)
            for snr_value in snr_values:
                detections_subset = detections[detections['SNR Value'] >= snr_value]
                relevant_column = detections_subset['Injected']
                collection[val][snr_value].append([np.sum(relevant_column), len(relevant_column) -
                                               np.sum(relevant_column)])


def contrast_curve(vals_and_finders, filepath_to_save):


    # Plotting Script From Planetfinder
    plt.plot(contrast_seps, contrast)
    plt.semilogy()
    plt.xlabel('Seperation (pixels)')
    plt.ylabel(r'$5\sigma$ contrast (Not Calibrated)')
    plt.title(r'$5\sigma$ Contrast at {0} $\mu$m (Not Calibrated)'.format(wavelength))
    plt.savefig('HD1160-KL20cc-{0}um-uncal.png'.format(wavelength))