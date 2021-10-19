from plotting_funcs import max_value_heatmap, paramvaluesfinder
import sys
import os

try:
    assert len(sys.argv) == 3  # the first argument in sys.argv is the name of the script ("maxheatmap.py")
except AssertionError:
    raise ValueError("Incorrect number of arguments. Specify two arguments: either 'annuli', 'subsections', "
                     "or 'movement'.")

# Using Log File to Get Param Values
annuli_vals = paramvaluesfinder('Annuli')
subsections_vals = paramvaluesfinder('Subsections')
movement_vals = paramvaluesfinder('Movement')

# Making Tuples to Pass Into max_value_heatmap
annuli = ('annuli', annuli_vals)
subsections = ('subsections', subsections_vals)
movement = ('movement', movement_vals)

# Output Filepaths
of1 = 'max_heatmap/ann-sbs-MAX.png'
of2 = 'max_heatmap/ann-mov-MAX.png'
of3 = 'max_heatmap/sbs-mov-MAX.png'

if not os.path.exists('max_heatmap'):
    try:
        os.mkdir('max_heatmap')
    except FileExistsError:  # if multiple scripts running at once (eg, if using analysis.sh)
        pass

params = [str.lower(param) for param in sys.argv[1:]]

if 'annuli' in params and 'subsections' in params:
    max_value_heatmap(annuli, subsections, of1)
elif 'annuli' in params and 'movement' in params:
    max_value_heatmap(annuli, movement, of2)
elif 'subsections' in params and 'movement' in params:
    max_value_heatmap(subsections, movement, of3)
else:
    raise ValueError(f'Sorry, script does not currently support params {params[0]} and {params[1]}')
