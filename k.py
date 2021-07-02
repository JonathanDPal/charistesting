import pandas as pd
import numpy as np

df = pd.read_csv('8Annuli_4Subsections_1Movement_NoneSpectrum_1Smooth_TrueHighpass__KL20_SNR-2.csv')

seps = [20, 40, 60]
pas = [0, 120, 240, 40, 160, 280, 80, 200, 320]
fwhm = 6.01924560185
masks = [[28, 52], [57, 37], [29, -30]]

def pasep_to_xy(PAs, seps):
	"""
	Takes lists of position angles and seperations and yields a numpy array with x-y coordinates for each combo, using
	the convention used in the table of values outputted by the planet detection software. The origin used in their
	convention is x=100, y=100 with respect to the bottom left corner of the image.
	"""
	radians = np.array(PAs) / 180 * np.pi
	locs = []
	for sep in seps:
		for rad in radians:
			x = -np.sin(rad) * sep
			y = np.cos(rad) * sep
			loc = [x, y]
			locs.append(loc)
	return locs

locs = pasep_to_xy(pas, seps)

def distance(xy1, xy2):
	return np.sqrt((xy1[0] - xy2[0]) ** 2 + (xy1[1] - xy2[1]) ** 2)

xpos = []
ypos = []
for x in df['x']:
	xpos.append(x)
for y in df['y']:
	ypos.append(y)
candidates = zip(xpos, ypos)

distances_from_fakes = []
distances_from_targets = []
for c in candidates:
	distances = []
	for loc in locs:
		distances.append(distance(c, loc))
	distances_from_fakes.append(np.min(distances))
	distances2 = []
	for mask in masks:
		distances2.append(distance(c, mask))
	distances_from_targets.append(np.min(distances2))

injected2 = []

	if d1 < fwhm:
		injected2.append("Science Target")
	elif d2 < fwhm:
		injected2.append(True)
	else:
		injected2.append(False)

df['Injected 2'] = injected2

print(df)
print(np.sum(df[df['Injected 2'] != "Science Target"]['Injected 2']))
