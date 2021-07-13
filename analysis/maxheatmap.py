from plotting_funcs import max_value_heatmap, paramvaluesfinder
import sys
import os

try:
    assert len(sys.argv) == 4
except AssertionError:
    raise ValueError("Incorrect number of arguments. Specify three arguments: first should be the name of the "
                     "object and the other two should be either 'annuli', 'subsections', or 'movement'.")

name = sys.argv[1]  # name of directory with data

# Using Log File to Get Param Values
annuli_vals = paramvaluesfinder(name, 'Annuli')
subsections_vals = paramvaluesfinder(name, 'Subsections')
movement_vals = paramvaluesfinder(name, 'Movement')

# Making Tuples to Pass Into max_value_heatmap
annuli = ('annuli', annuli_vals)
subsections = ('subsections', subsections_vals)
movement = ('movement', movement_vals)

# Output Filepaths
of1 = 'max_heatmap/ann-sbs.png'
of2 = 'max_heatmap/ann-mov.png'
of3 = 'max_heatmap/sbs-mov.png'

if not os.path.exists('max_heatmap'):
    os.mkdir('max_heatmap')

params = [str.lower(param) for param in sys.argv[2:4]]

if 'annuli' in params and 'subsections' in params:
        max_value_heatmap(annuli, subsections, name, of1)
elif 'annuli' in params and 'movement' in params:
    max_value_heatmap(annuli, movement, name, of2)
elif 'subsections' in params and 'movement' in params:
    max_value_heatmap(subsections, movement, name, of3)
else:
    raise ValueError(f'Sorry, script does not currently support params {params[0]} and {params[1]}')
