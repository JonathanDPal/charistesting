1. Get approximate FWHMs for all the different filters
and add those in to the FWHMIOWA_calculator function
in parameter_test_infrastructure.py. The actual
implementation has already been outlined in comments
in that function.

2. Add another script in organzational_scripts that gives
you optimal planet injection sites based on how many planets
you want to inject, the annuli/subsections values you want to
use, and the locations of the satellite spots/science targets
in the observations you are going to use.

3. Add some stuff so that you have more flexibility on
fake_seps.

4. Add some stuff so that if fake planets would otherwise
be injected on annuli/subsection boundaries, they get
tweaked slighly automatically.

5. Test the implementation for the TODO items which have
been marked as complete but not yet tested.