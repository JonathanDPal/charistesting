from plotting_funcs import max_value_heatmap

annuli = ('annuli', [4, 6, 8, 10, 12])
subsections = ('subsections', [2, 4, 6])
movement = ('movement', [0, 1, 2])

ds = '../detections/'
of = 'max_heatmap/ann-sbs'

max_value_heatmap(annuli, subsections, ds, of)


