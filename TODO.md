1. Get approximate FWHMs for all the different filters
and add those in to the FWHMIOWA_calculator function
in parameter_test_infrastructure.py. The actual
implementation has already been outlined in comments
in that function.

2. Add some stuff so that if fake planets would otherwise
be injected on annuli/subsection boundaries, they get
tweaked slighly automatically.

3. Test the implementation for the TODO items which have
been marked as complete but not yet tested.