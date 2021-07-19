from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from pyklip.fakes import retrieve_planet_flux
from copy import copy
from parameter_test_infrastructure import make_dn_per_contrast, FWHMIOWA_calculator
from pyklip.instruments.CHARIS import CHARISData
import pandas as pd

orig_filepath = 'HD1160_cubes/CRSA00015235_cube.fits'
dataset = make_dn_per_contrast(CHARISData(orig_filepath))
dn_per_contrast = dataset.dn_per_contrast
dataset = None

filepath = 'HD1160/klipped_cubes_Wfakes' \
            '/HD1160_withfakes_4Annuli_2Subsections_0Movement_NoneSpectrum_0Smooth_TrueHighpass_-KL20-speccube.fits'
with fits.open(filepath) as hdulist:
    cube = copy(hdulist[1].data)
    dataset_center = [hdulist[1].header['PSFCENTX'], hdulist[1].header['PSFCENTY']]
    dataset_fwhm, dataset_iwa, dataset_owa = FWHMIOWA_calculator(hdulist)
    output_wcs = WCS(hdulist[0].header, naxis=[1, 2])

wavelength_index = 10
frame = cube[wavelength_index] / dn_per_contrast[wavelength_index]

mask_xy = [144, 80]

x_pos = mask_xy[0]
y_pos = mask_xy[1]
ydat, xdat = np.indices(frame.shape)
distance_from_planet = np.sqrt((xdat - x_pos) ** 2 + (ydat - y_pos) ** 2)
frame[np.where(distance_from_planet <= 2 * dataset_fwhm)] = np.nan

fake_fluxes = [5e-4, 5e-5, 5e-6, 1e-4, 1e-5, 1e-6]  # List of Float(s)
fake_seps = [20, 40, 60]  # List of Integer(s) and/or Float(s)
fake_PAs = [19, 79, 139, 199, 259, 319]  # List of Integer(s) and/or Float(s)

retrieved_fluxes = []
for sep in fake_seps:
    for pa in fake_PAs:
        fake_flux = retrieve_planet_flux(frame, dataset_center, output_wcs, sep, pa, guessfwhm=dataset_fwhm,
                                         guesspeak=1e-6, refinefit=True, searchrad=12)
        retrieved_fluxes.append(fake_flux)

at20 = [fake_fluxes[0], fake_fluxes[3]] * 3
at20.append(float(np.mean(at20)))
at40 = [fake_fluxes[1], fake_fluxes[4]] * 3
at40.append(float(np.mean(at40)))
at60 = [fake_fluxes[2], fake_fluxes[5]] * 3
at60.append(float(np.mean(at60)))

assert len(retrieved_fluxes) == 18

act20 = retrieved_fluxes[:6] + [np.mean(retrieved_fluxes[:6])]
act40 = retrieved_fluxes[6:12] + [np.mean(retrieved_fluxes[6:12])]
act60 = retrieved_fluxes[12:] + [np.mean(retrieved_fluxes[12:])]

df = pd.DataFrame(index=([str(num+1) for num in range(6)] + ['MEAN']))

df['Injected (20 sep)'] = at20
df['Actual (20 sep)'] = act20
df['Injected (40 sep)'] = at40
df['Actual (40 sep)'] = act40
df['Injected (60 sep)'] = at60
df['Actual (60 sep)'] = act60

df.to_csv('fluxes_refinefitTrue.csv')
