import pylab as P
import pandas as pd
import numpy as np

snr_values = np.arange(start=3, stop=10, step=0.5)
data = dict()
for snr_value in snr_values:
    data[snr_value] = list()

file_finder = ''
fileset = glob(file_finder)
for file in fileset:
    detections = pd.read_csv('')
    for snr_value in snr_values:
        detections_subset = detections[detections['SNR Value'] >= snr_value]
        relevant_column = detections_subset['Injected']
        data[snr_value].append([np.sum(relevant_column), len(relevant_column) - np.sum(
            relevant_column)])
