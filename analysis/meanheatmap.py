import sys
from plotting_funcs import mean_value_heatmap, paramvaluesfinder
import os

try:
    assert len(sys.argv) == 4
except AssertionError:
    raise ValueError("Incorrect number of arguments. First argument should be the name of the object. Second and third "
                     "arguments should be either 'annuli', 'subsections', or 'movement'.")

name = sys.argv[1]  # name of directory with all of the data in it

# Using Log File to Get Values
annuli_vals = paramvaluesfinder(name, 'Annuli')
subsections_vals = paramvaluesfinder(name, 'Subsections')
movement_vals = paramvaluesfinder(name, 'Movement')
ni = paramvaluesfinder(name, 'ni')  # number of injected planets

# Putting Into Tuples For mean_value_heatmap
annuli = ('annuli', annuli_vals)
subsections = ('subsections', subsections_vals)
movement = ('movement', movement_vals)

# Output Filepaths
of1 = 'max_heatmap/ann-sbs'
of2 = 'max_heatmap/ann-mov'
of3 = 'max_heatmap/sbs-mov'

if not os.path.exists('max_heatmap'):
    os.mkdir('max_heatmap')

params = [str.lower(param) for param in sys.argv[2:4]]

if 'annuli' in params and 'subsections' in params:
    mean_value_heatmap(annuli, subsections, ni, name, of1)
elif 'annuli' in params and 'movement' in params:
    mean_value_heatmap(annuli, movement, ni, name, of2)
elif 'subsections' in params and 'movement' in params:
    mean_value_heatmap(subsections, movement, ni, name, of3)
else:
    raise ValueError(f'Sorry, script does not currently support params {params[0]} and {params[1]}')
