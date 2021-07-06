import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from glob import glob
import os

def valuefinder(filename, param):
    param = str.lower(param)
    paramlengths = {'annuli': 6, 'subsections': 11, 'movement': 8, 'spectrum': 8, 'smooth': 6, 'highpass': 8, 'KL': 4}
    paramlength = paramlengths[param][0]
    startingindex = None # will be defined soon
    for i in range(len(filename)):
        if str.lower(filename[i:i+paramlength]) == param:
            startingindex = i - 1

    valuelength = 0
    while startingindex >= 0 and filename[startingindex] != '_':
        startingindex -= 1
        valuelength += 1
    value = filename[startingindex:startingindex+valuelength]

    return value


def roc_generator(snr_values, vals_and_finders, filepath_to_save):
    """
    Makes an ROC curve.
    ---
    Args:
        snr_values (list): SNR Values to plot
        vals_and_finders (dict): Values to check out are keys and file finders are the values.
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
                detections_subset = detections[detections['SNR Value'] >= snr_value and detections['Injected'] !=
                                               "Science Target"]
                detections['Injected'] = [True if x == "True" else False for x in detections['Injected']]
                inj = detections_subset['Injected']
                collection[val][snr_value].append([np.sum(inj), len(inj) - np.sum(inj)])

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


def max_value_heatmap(param1, param2, directory_switcher, filepath_to_save, file_finder='*.csv'):
    """
    Just shows the maximum SNR value found (i.e. SNR of HD1160)
    Args:
        param1: Should be tuple of form (str: name of parameter, list: values used for parameter)
        param2: Should be tuple of form (str: name of parameter, list: values used for parameter)
        directory_switcher: string that can get passed in to change working directory to "detections"
        filepath_to_save: string
        file_finder: str -- passed into glob to get all relevant (CSV) files.
    """
    os.chdir(directory_switcher)
    fileset = glob(file_finder)
    full_data = {str(a): {str(m):[] for m in param2[1]} for a in param1[1]}
    for file in fileset:
        df = pd.read_csv(file)
        snr = df['SNR Value'].max()
        p1 = valuefinder(file, param1[0])
        p2 = valuefinder(file, param2[0])
        full_data[p1][p2].append(snr)
    for p1 in param1[1]:
        for p2 in param2[1]:
            full_data[p1][p2] = np.mean(full_data[p1][p2])
    plot_snr = []
    for p2 in param2[1]:
        plot_snr.append([full_data[p1][p2] for p1 in param1[1]])
    data_to_plot = pd.DataFrame(plot_snr, index=param2, columns=param1[1])
    sns.heatmap(data_to_plot, annot=True, linewidths=0.2, fmt='d')
    plt.savefig(filepath_to_save)


def mean_value_heatmap(param1, param2, num_injections, directory_switcher, filepath_to_save, file_finder='*.csv'):
    """
    Shows mean SNR of injected planets.

    Args:
        param1: Should be tuple of form (str: name of parameter, list: values used for parameter)
        param2: Should be tuple of form (str: name of parameter, list: values used for parameter)
        num_injections: Number of fake planets injected.
        directory_switcher: string that can get passed in to change working directory to "detections"
        filepath_to_save: string
        file_finder: str -- passed into glob to get all relevant (CSV) files.
    """
    os.chdir(directory_switcher)
    fileset = glob(file_finder)
    full_data = {str(a): {str(m):[] for m in param2[1]} for a in param1[1]}
    for file in fileset:
        df = pd.read_csv(file)
        injected = df[df["Injected"] == "True"]
        snr = injected['SNR Value'].mean()
        p1 = valuefinder(file, param1[0])
        p2 = valuefinder(file, param2[0])
        full_data[p1][p2].append(snr)
    for p1 in param1[1]:
        for p2 in param2[1]:
            if not len(full_data[p1][p2]) == num_injections:
                missing = [1] * (len(full_data[p1][p2]) - num_injections)
                full_data[p1][p2].append(missing)
            full_data[p1][p2] = np.mean(full_data[p1][p2])
    plot_snr = []
    for p2 in param2[1]:
        plot_snr.append([full_data[p1][p2] for p1 in param1[1]])
    data_to_plot = pd.DataFrame(plot_snr, index=param2, columns=param1[1])
    sns.heatmap(data_to_plot, annot=True, linewidths=0.2, fmt='d')
    plt.savefig(filepath_to_save)





