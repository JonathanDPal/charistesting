import numpy as np
import pandas as pd

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

def distance(xy1, xy2):
	return np.sqrt((xy1[0] - xy2[0]) ** 2 + (xy1[1] - xy2[1]) ** 2)

for i in range(40):
	pas = [0, 60, 120, 180, 240, 300]
	pas = np.array(pas) + i
	seps = [20, 40, 60]

	locs = pasep_to_xy(pas, seps)

	hd1160 = pd.read_excel('satellite_spots.xlsx', sheet_name='HD1160')
	hr8799 = pd.read_excel('satellite_spots.xlsx', sheet_name='HR8799')
	masks = [[28, 52], [57, 37], [29, -30], [44, -20]]

	xpos = []
	ypos = []
	for x in hd1160['KLIP x-pos']:
		xpos.append(x - 100)
	for y in hd1160['KLIP y-pos']:
		ypos.append(y - 100)
	for x in hr8799['KLIP x-pos']:
		xpos.append(x - 100)
	for y in hr8799['KLIP y-pos']:
		ypos.append(y - 100)
	spots = zip(xpos, ypos)


	distances = []
	for loc in locs:
		for spot in spots:
			distances.append(distance(loc, spot))
		for mask in masks:
			distances.append(distance(loc, mask))

	print(i, " : ", np.min(distances))

def FWHMIOWA_calculator(speccubefile, filtname=None):
	"""
	Finds FWHM, IWA, and OWA for a opened CHARIS data cube. Thanks to Dr. Tobin for this.
	(https://docs.google.com/document/d/1S1Oo9QweKwnOfmv6fu28bb75lYeQXGzn/edit)
	"""
	wavelengths = {'j': 1200e-9, 'h': 1550e-9, 'k': 2346e-9, 'broadband': 1550e-9}
	if filtname is None:
		wavelength = wavelengths[str.lower(speccubefile[1].header['FILTNAME'])]
	else:
		wavelength = wavelengths[str.lower(filtname)]
	D = 8
	lenslet_scale = 0.0162
	field_radius = 1.035
	FWHM = 2 * 1.22 * wavelength / D * 206265 / lenslet_scale
	IWA = 5
	OWA = (field_radius / lenslet_scale) - FWHM

	return FWHM, IWA, OWA

print(FWHMIOWA_calculator(speccubefile=None, filtname='Broadbandr')[0] * 2)
print(FWHMIOWA_calculator(speccubefile=None, filtname='K')[0] * 2)