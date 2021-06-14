from test_parameter_test_infrastructure import *

# Measuring Contrast, Generating Contrast Curve on Data Without Fakes &
# Saving the Contrast Data to a CSV File
no_fakes_wln_um, no_fakes_calib_cube, no_fakes_dataset_center, no_fakes_dataset_iwa, \
no_fakes_dataset_owa, no_fakes_dataset_fwhm, no_fakes_output_wcs = ss_calibration(False)

get_contrast(graph_output_filepath=filepath_uncalibrated_curve,
						data_output_filepath=filepath_uncal_contrast_data,
						wln_um=no_fakes_wln_um, calib_cube=no_fakes_calib_cube,
						dataset_center=no_fakes_dataset_center,
						dataset_iwa=no_fakes_dataset_iwa, dataset_owa=no_fakes_dataset_owa,
						dataset_fwhm=no_fakes_dataset_fwhm, output_wcs=no_fakes_output_wcs)

# Measuring Contrast, Generating Contrast Curve on Data With Fakes &
# Saving the Contrast Data to a CSV File
with_fakes_wln_um, with_fakes_calib_cube, with_fakes_dataset_center, \
with_fakes_dataset_iwa, with_fakes_dataset_owa, with_fakes_dataset_fwhm, \
with_fakes_output_wcs = ss_calibration(True)

get_contrast(graph_output_filepath=filepath_calibrated_curve,
						data_output_filepath=filepath_cal_contrast_data,
						wln_um=with_fakes_wln_um, calib_cube=with_fakes_calib_cube,
						dataset_center=with_fakes_dataset_center,
						dataset_iwa=with_fakes_dataset_iwa,
						dataset_owa=with_fakes_dataset_owa,
						dataset_fwhm=with_fakes_dataset_fwhm,
						output_wcs=with_fakes_output_wcs, contains_fakes=True,
						injected_fluxes=self.fake_fluxes, injected_seps=self.fake_seps,
						injected_PAs=self.fake_PAs)