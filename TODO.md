1. Write the scoring algorithm and implement it in
a Python script.

2. Get approximate FWHMs for all the different filters
and add those in to the FWHMIOWA_calculator function
in parameter_test_infrastructure.py. The actual
implementation here has already been outlined in comments
in that function.

3. Finish writing the injection_locations.py script

4. In addition to 2, add another script in 
organzational_scripts that gives you optimal planet 
injection sites based on how many planets you want to 
inject, the annuli/subsections values you want to use,
and the locations of the satellite spots/science targets
in the observations you are going to use.

5. Add some stuff so that you have more flexibility on
fake_seps.

6. Add some stuff so that if fake planets would otherwise
be injected on annuli/subsection boundaries, they get
tweaked slighly automatically.

7. Add a _verbose_ attribute for TestDataset so that log
files can be super decluttered if people want that.
