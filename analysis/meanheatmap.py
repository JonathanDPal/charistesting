import sys
from plotting_funcs import mean_value_heatmap, paramvaluesfinder
import os

try:
    assert len(sys.argv) == 3  # the first arg in sys.argv is the name of the script ("meanheatmap.py")
except AssertionError:
    raise ValueError("Incorrect number of keyword arguments. Specify two arguments which should be either 'annuli', "
                     "'subsections', or 'movement'.")

# Using Log File to Get Values
annuli_vals = paramvaluesfinder('Annuli')
subsections_vals = paramvaluesfinder('Subsections')
movement_vals = paramvaluesfinder('Movement')
ni = paramvaluesfinder('ni')  # number of injected planets

# Putting Into Tuples For mean_value_heatmap
annuli = ('annuli', annuli_vals)
subsections = ('subsections', subsections_vals)
movement = ('movement', movement_vals)

# Output Filepaths
of1 = 'mean_heatmap/ann-sbs'
of2 = 'mean_heatmap/ann-mov'
of3 = 'mean_heatmap/sbs-mov'

if not os.path.exists('mean_heatmap'):
    try:
        os.mkdir('mean_heatmap')
    except FileExistsError:  # if multiple scripts running at once (eg, if using analysis.sh)
        pass

params = [str.lower(param) for param in sys.argv[1:]]

if 'annuli' in params and 'subsections' in params:
    mean_value_heatmap(annuli, subsections, ni, of1)
elif 'annuli' in params and 'movement' in params:
    mean_value_heatmap(annuli, movement, ni, of2)
elif 'subsections' in params and 'movement' in params:
    mean_value_heatmap(subsections, movement, ni, of3)
else:
    raise ValueError(f'Sorry, script does not currently support params {params[0]} and {params[1]}')
