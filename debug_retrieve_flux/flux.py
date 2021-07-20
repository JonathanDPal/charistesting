from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from pyklip.fakes import retrieve_planet_flux, retrieve_planet
from copy import copy
from parameter_test_infrastructure import make_dn_per_contrast, FWHMIOWA_calculator, pasep_to_xy
from pyklip.instruments.CHARIS import CHARISData
import pandas as pd

orig_filepath = 'HD1160_cubes/CRSA00015235_cube.fits'
dataset = make_dn_per_contrast(CHARISData(orig_filepath))
dn_per_contrast = dataset.dn_per_contrast
dataset = None

filepath = 'HD1160/klipped_cubes_Wfakes' \
            '/HD1160_withfakes_4Annuli_2Subsections_0.0Movement_NoneSpectrum_0.0Smooth_TrueHighpass_-KL20-speccube.fits'
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
ret_fwhms = []
ret_xfits = []
ret_yfits = []
for sep in fake_seps:
    for pa in fake_PAs:
        fake_flux, xfit, yfit, ret_fwhm = retrieve_planet(frame, dataset_center, output_wcs, sep, pa,
                                                           guessfwhm=dataset_fwhm, guesspeak=1e-6, refinefit=True,
                                                           searchrad=12)
        retrieved_fluxes.append(fake_flux)
        ret_fwhms.append(ret_fwhm)
        ret_xfits.append(xfit - dataset_center[0])
        ret_yfits.append(yfit - dataset_center[0])

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

fwhms20 = ret_fwhms[:6] + [np.mean(ret_fwhms[:6])]
fwhms40 = ret_fwhms[6:12] + [np.mean(ret_fwhms[6:12])]
fwhms60 = ret_fwhms[12:] + [np.mean(ret_fwhms[12:])]

xfits20 = ret_xfits[:6] + [np.nan]
xfits40 = ret_xfits[6:12] + [np.nan]
xfits60 = ret_xfits[12:] + [np.nan]

yfits20 = ret_yfits[:6] + [np.nan]
yfits40 = ret_yfits[6:12] + [np.nan]
yfits60 = ret_yfits[12:] + [np.nan]

locs = pasep_to_xy(fake_PAs, fake_seps)

xlocs20 = [loc[0] for loc in locs[:6]] + [np.nan]
xlocs40 = [loc[0] for loc in locs[6:12]] + [np.nan]
xlocs60 = [loc[0] for loc in locs[12:]] + [np.nan]

ylocs20 = [loc[1] for loc in locs[:6]] + [np.nan]
ylocs40 = [loc[1] for loc in locs[6:12]] + [np.nan]
ylocs60 = [loc[1] for loc in locs[12:]] + [np.nan]

df = pd.DataFrame(index=([str(num+1) for num in range(6)] + ['MEAN']))

df['Injected Flux (20 sep)'] = at20
df['Retrieved Flux (20 sep)'] = act20
df['Injected X (20 sep)'] = xlocs20
df['Retrieved X (20 sep)'] = xfits20
df['Injected Y (20 sep)'] = ylocs20
df['Retrieved Y (20 sep)'] = yfits20
df['Retrieved FWHM (20 sep)'] = fwhms20

df['Injected Flux (40 sep)'] = at40
df['Retrieved Flux (40 sep)'] = act40
df['Injected X (40 sep)'] = xlocs40
df['Retrieved X (40 sep)'] = xfits40
df['Injected Y (40 sep)'] = ylocs40
df['Retrieved Y (40 sep)'] = yfits40
df['Retrieved FWHM (40 sep)'] = fwhms40

df['Injected Flux (60 sep)'] = at60
df['Retrieved Flux (60 sep)'] = act60
df['Injected X (60 sep)'] = xlocs60
df['Retrieved X (60 sep)'] = xfits60
df['Injected Y (60 sep)'] = ylocs60
df['Retrieved Y (60 sep)'] = yfits60
df['Retrieved FWHM (60 sep)'] = fwhms60

df.to_excel('retrieve_planet_data.xlsx')
