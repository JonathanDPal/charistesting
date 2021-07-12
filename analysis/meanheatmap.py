from plotting_funcs import mean_value_heatmap

annuli = ('annuli', [4, 6, 8, 10, 12])
subsections = ('subsections', [2, 4, 6])
movement = ('movement', [0, 1, 2])

ds = '../detections'
of = 'mean_heatmap/ann-sbs'

ni = 18

mean_value_heatmap(annuli, subsections, ni, ds, of)
