import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import pandas as pd
import numpy as np
from glob import glob
import os


def paramvaluesfinder(object_name, param):
    paramline = None
    fluxes = None
    pas = None
    with open('/home/jpal/data0/jpal/parameter_sampling/' + object_name + '/log.txt') as logfile:
        for line in logfile:
            if str.lower(param) != 'ni':
                if str.lower(param) in str.lower(line):
                    paramline = line
                    break
            else:
                if str.lower('Fake Fluxes') in str.lower(line):
                    fluxes = line
                elif str.lower('Fake PAs') in str.lower(line):
                    pas = line
    if paramline is not None:
        paramline = paramline.replace(' ', '')
        for i in range(len(paramline)):
            if paramline[i] == '[':
                starting_index = i + 1
            elif paramline[i] == ']':
                final_index = i
        if str.lower(param) in ['annuli', 'subsections', 'numbasis', 'movement', 'corr_smooth']:
            vals = [int(val) for val in paramline[starting_index: final_index].split(',')]
        elif str.lower(param) in ['none at the moment']:
            vals = [float(val) for val in paramline[starting_index: final_index].split(',')]
        elif str.lower(param) in ['highpass']:
            vals = list()
            for val in paramline[starting_index: final_index].split(','):
                if str.lower(val) == 'true':
                    vals.append(True)
                elif str.lower(val) == 'false':
                    vals.append(False)
                else:
                    vals.append(float(val))
        else:
            raise ValueError(f"Sorry, this function does not currently support value finding for param {param}.")
    elif fluxes is not None and pas is not None:
        flux_vals = list()
        pa_vals = list()
        fluxes = fluxes.replace(' ', '')
        pas = pas.replace(' ', '')
        for p in [fluxes, pas]:
            for i in range(len(p)):
                if p[i] == '[':
                    starting_index = i + 1
                elif p[i] == ']':
                    final_index = i
            for val in p[starting_index: final_index].split(','):
                if p == fluxes:
                    flux_vals.append(float(val))
                elif p == pas:
                    pa_vals.append(float(val))
        vals = len(flux_vals) * len(pa_vals)

    return vals


def valuefinder(filename, param):
    param = str.lower(param)
    paramlengths = {'annuli': 6, 'subsections': 11, 'movement': 8, 'spectrum': 8, 'smooth': 6, 'highpass': 8, 'kl': 2}
    paramlength = paramlengths[param]
    startingindex = None  # will be defined soon
    for i in range(len(filename)):
        if str.lower(filename[i: i + paramlength]) == param:
            if param == 'kl':
                startingindex = i + paramlength
            else:
                startingindex = i - 1

    if param != 'kl':
        valuelength = 0
        while startingindex >= 0 and filename[startingindex] != '_':
            startingindex -= 1
            valuelength += 1
        value = filename[startingindex + 1: startingindex + 1 + valuelength]
    else:
        value = filename[startingindex: startingindex + 2]

    return value


def roc_generator(snr_values, param1, num_injections, object_name, filepath_to_save):
    """
    Makes an ROC curve.
    ---
    Args:
        snr_values (list): SNR Values to plot
        param1: Should be tuple of form (str: name of parameter, list: values used for parameter)
        num_injections (int): Number of planets that got injected
        object_name (str)
        filepath_to_save (str): Where to save graph.
    """
    collection = {str(val): {str(snr): list() for snr in snr_values} for val in param1[1]}
    originalwd = os.getcwd()
    os.chdir('/home/jpal/data0/jpal/parameter_sampling/' + object_name + '/detections/')
    filelist = glob('*.csv')
    for file in filelist:  # i used only to get different colors/markers on graph
        detections = pd.read_csv(file)
        val = valuefinder(file, param1[0])
        for snr in snr_values:
            detections_subset = detections[detections['SNR Value'] >= snr]
            detections_subset = detections_subset[detections_subset['Injected'] != "Science Target"]
            inj = []
            for m in detections_subset['Injected']:
                if m is True or m == "True":
                    inj.append(True)
                elif m is False or m == "False":
                    inj.append(False)
            collection[val][str(snr)].append([np.sum(inj), len(inj) - np.sum(inj)])
    for i, val in enumerate(param1[1]):
        val = str(val)
        k = collection[val]
        x = []
        y = []
        for snr in snr_values:
            A = k[str(snr)]
            tp = np.sum([a[0] for a in A])
            fp = np.sum([a[1] for a in A])
            y.append(tp / (num_injections * len(A)))
            x.append(fp / (num_injections * len(A)))
        markers = ['.', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', 'P', '*', 'h', 'H', '+',
                   'x', 'X', 'D', 'd', '|', '_']
        colors = list(mcolors.BASE_COLORS) + list(mcolors.TABLEAU_COLORS)
        plt.plot(x, y, label=val, marker=markers[i], color=colors[i])

    plt.xlabel('False Positives')
    plt.ylabel('True Positives')
    plt.legend(loc='lower right')
    os.chdir(originalwd)
    plt.savefig(filepath_to_save)


def max_value_heatmap(param1, param2, object_name, filepath_to_save, file_finder='*.csv'):
    """
    Just shows the maximum SNR value found (i.e. SNR of HD1160)
    Args:
        param1: Should be tuple of form (str: name of parameter, list: values used for parameter)
        param2: Should be tuple of form (str: name of parameter, list: values used for parameter)
        object_name
        filepath_to_save: string
        file_finder: str -- passed into glob to get all relevant (CSV) files.
    """
    originalwd = os.getcwd()
    os.chdir('/home/jpal/data0/jpal/parameter_sampling/' + object_name + '/detections/')
    fileset = glob(file_finder)
    full_data = {str(a): {str(m): list() for m in param2[1]} for a in param1[1]}
    for file in fileset:
        df = pd.read_csv(file)
        snr = df['SNR Value'].max()
        p1 = valuefinder(file, param1[0])
        p2 = valuefinder(file, param2[0])
        full_data[p1][p2].append(snr)
    for p1 in param1[1]:
        for p2 in param2[1]:
            p1 = str(p1)
            p2 = str(p2)
            full_data[p1][p2] = np.mean(full_data[p1][p2])
            full_data[p1][p2] = int(round(full_data[p1][p2]))
    plot_snr = []
    for p2 in param2[1]:
        plot_snr.append([full_data[str(p1)][str(p2)] for p1 in param1[1]])
    data_to_plot = pd.DataFrame(plot_snr, index=param2[1], columns=param1[1])
    sns.heatmap(data_to_plot, annot=True, linewidths=0.2, fmt='d')
    os.chdir(originalwd)
    plt.savefig(filepath_to_save)


def mean_value_heatmap(param1, param2, num_injections, object_name, filepath_to_save, file_finder='*.csv'):
    """
    Shows mean SNR of injected planets.
    Args:
        param1: Should be tuple of form (str: name of parameter, list: values used for parameter)
        param2: Should be tuple of form (str: name of parameter, list: values used for parameter)
        num_injections: Number of fake planets injected.
        object_name
        filepath_to_save: string
        file_finder: str -- passed into glob to get all relevant (CSV) files.
    """
    originalwd = os.getcwd()
    os.chdir('/home/jpal/data0/jpal/parameter_sampling/' + object_name + '/detections/')
    fileset = glob(file_finder)
    full_data = {str(a): {str(m): list() for m in param2[1]} for a in param1[1]}
    for file in fileset:
        df = pd.read_csv(file)
        injected = df[df["Injected"] == "True"]
        snrs = list(injected['SNR Value'])
        missing = [0] * (len(injected['SNR Value']) - num_injections)  # heavily punishing non-detections
        snr_avg = np.mean(snrs + missing)
        p1 = valuefinder(file, param1[0])
        p2 = valuefinder(file, param2[0])
        try:
            full_data[p1][p2].append(snr_avg)
        except KeyError:
            print(full_data)
            print(p1, p2)
            raise KeyError
    for p1 in param1[1]:
        for p2 in param2[1]:
            p01 = str(p1)
            p02 = str(p2)
            full_data[p01][p02] = np.nanmean(full_data[p01][p02])
            full_data[p01][p02] = round(full_data[p01][p02], 2)
    plot_snr = []
    for p2 in param2[1]:
        plot_snr.append([full_data[str(p1)][str(p2)] for p1 in param1[1]])
    data_to_plot = pd.DataFrame(plot_snr, index=param2[1], columns=param1[1])
    sns.heatmap(data_to_plot, annot=True, linewidths=0.2, fmt='.2f')
    os.chdir(originalwd)
    plt.savefig(filepath_to_save)
