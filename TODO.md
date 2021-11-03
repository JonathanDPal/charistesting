1. Get approximate FWHMs for all the different filters
and add those in to the FWHMIOWA_calculator function
in parameter_test_infrastructure.py. The actual
implementation here has already been outlined in comments
in that function.

2. Finish writing the injection_locations.py script

3. In addition to 2, add another script in 
organzational_scripts that gives you optimal planet 
injection sites based on how many planets you want to 
inject, the annuli/subsections values you want to use,
and the locations of the satellite spots/science targets
in the observations you are going to use.

4. Add some stuff so that you have more flexibility on
fake_seps.

5. Add some stuff so that if fake planets would otherwise
be injected on annuli/subsection boundaries, they get
tweaked slighly automatically.

6. Add a _verbose_ attribute for TestDataset so that log
files can be super decluttered if people want that.

7. Add an argument where people can disable the creation
of a log file (in case just doing some quick fill in
and don't want old log file overwritten)

8. Fix paramvaluesfinder in plotting_funcs.py such that it
can accurately identify the number of planets injected if
the flux/sep lists are lists of lists