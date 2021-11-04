from plotting_funcs import specific_target_heatmap, paramvaluesfinder
import sys
import os

# Using Log File to Get Param Values
annuli_vals = paramvaluesfinder('Annuli')
subsections_vals = paramvaluesfinder('Subsections')
movement_vals = paramvaluesfinder('Movement')

# Making Tuples to Pass Into max_value_heatmap
annuli = ('annuli', annuli_vals)
subsections = ('subsections', subsections_vals)
movement = ('movement', movement_vals)

# Locations of Planets Within Images
loc_d = [129, 70]
loc_c = [128, 152]
loc_b = [157, 137]

# Output Filepaths
ofd0 = 'target_heatmaps/ann-sbs-D.png'
ofd1 = 'target_heatmaps/ann-mov-D.png'
ofd2 = 'target_heatmaps/sbs-mov-D.png'
ofc0 = 'target_heatmaps/ann-sbs-C.png'
ofc1 = 'target_heatmaps/ann-mov-C.png'
ofc2 = 'target_heatmaps/sbs-mov-C.png'
ofb0 = 'target_heatmaps/ann-sbs-B.png'
ofb1 = 'target_heatmaps/ann-mov-B.png'
ofb2 = 'target_heatmaps/sbs-mov-B.png'

if not os.path.exists('target_heatmaps'):
    try:
        os.mkdir('target_heatmaps')
    except FileExistsError:  # if multiple scripts running at once (eg, if using analysis.sh)
        pass

params = [str.lower(param) for param in sys.argv[1:]]

if ('annuli' in params and 'subsections' in params) or 'all' in params:
    specific_target_heatmap(annuli, subsections, ofd0, loc_d, target_name='HR8799d')
    specific_target_heatmap(annuli, subsections, ofc0, loc_c, target_name='HR8799c')
    specific_target_heatmap(annuli, subsections, ofb0, loc_b, target_name='HR8799b')
if ('annuli' in params and 'movement' in params) or 'all' in params:
    specific_target_heatmap(annuli, movement, ofd1, loc_d, target_name='HR8799d')
    specific_target_heatmap(annuli, movement, ofc1, loc_c, target_name='HR8799c')
    specific_target_heatmap(annuli, movement, ofb0, loc_b, target_name='HR8799b')
if ('subsections' in params and 'movement' in params) or 'all' in params:
    specific_target_heatmap(subsections, movement, ofd2, loc_d, target_name='HR8799d')
    specific_target_heatmap(subsections, movement, ofc2, loc_c, target_name='HR8799c')
    specific_target_heatmap(subsections, movement, ofb2, loc_b, target_name='HR8799b')
