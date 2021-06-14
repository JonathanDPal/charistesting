import pytest
from parameter_test_infrastructure import *

speccubefile = 'test_files/HD1160-KL20-speccube.fits'

def test_FWHMIOWA_calculator_broadband(speccubefile):
    with fits.open(speccubefile) as hdulist:
        assert FWHMIOWA_calculator(hdulist) == fwhm, iwa, owa

def test_calibrate_ss_contrast(speccubefile):
    with fits.open(speccubefile) as hdulist:
        assert calibrate_ss_contrast(hdulist) == wln_um, spot_to_star

fake_fluxes = [5.0e-3, 1.0e-4]
fake_seps = [15, 30]
annuli = [5, 10]
subsections = [2, 4]
movement = [1, 2]
numbasis = [3, 5]
spectrum = ['methane', None]
fake_PAs =
HD1160 = TestDataset(fileset='test_files/hd1160/*.fits', object_name='HD_11_60', mask_xy=[144,80],
                     fake_fluxes=fake_fluxes, fake_seps=fake_seps, annuli=annuli,
                     subsections=subsections, movement=movement, numbasis=numbasis,
                     spectrum=spectrum, mode='ADI+SDI', fake_fwhm=3.5, fake_PAs=fake_PAs)

def test_trial_list_build():
    pass

def test_inject_fakes():
    pass

def test_run_KLIP():
    pass

def test_klip_parameter_string_creation():
    pass

def test_klipped_filepath_string_creation():
    pass

def test_get_contrast():
    pass

def test_detect_planets():
    pass

def test_categorize_detections():
    pass

