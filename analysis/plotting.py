import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from glob import glob
import os

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

    for val in vals_and_finders.keys():
        k = collection[val]
        x = []
        y = []
        for snr in snr_values:
            tp = np.sum([a[0] for a in k[snr]])
            y.append(tp)
            fp = np.sum([a[1] for a in k[snr]])
            x.append(fp)
        plt.plot(x,y, label=val)

    plt.xlabel('False Positives')
    plt.ylabel('True Positives')
    plt.legend(loc='upper right')
    plt.savefig(filepath_to_save)


def contrast_curve(vals_and_files, filepath_to_save):
    for val in vals_and_files.keys():
        data = pd.read_csv(vals_and_files[val])
        data.plot(label=val)

    plt.semilogy()
    plt.xlabel('Seperation (pixels)')
    plt.ylabel('Contrast')
    plt.legend(loc='upper right')
    plt.savefig(filepath_to_save)


def heatmap(annuli, movement, directory_switcher, filepath_to_save, file_finder='*.csv'):
    """
    For now, just shows the maximum SNR value found
    """
    os.chdir(directory_switcher)
    fileset = glob(file_finder)
    full_data = {a: {m:[] for m in movement} for a in annuli}
    for file in fileset:
        df = pd.read_csv(file)
        snr = df['SNR Value'].max()
        ani = int(file[0])
        mov = int(file[21])
        full_data[ani][mov].append(snr)
    for ani in annuli:
        for mov in movement:
            full_data[ani][mov] = np.mean(full_data[ani][mov])
    plot_snr = []
    for mov in movement:
        plot_snr.append([full_data[ani][mov] for ani in annuli])
    data_to_plot = pd.DataFrame(plot_snr, index=movement, columns=annuli)
    sns.heatmap(data_to_plot, annot=True, linewidths=0.2, fmt='d')
    plt.savefig(filepath_to_save)