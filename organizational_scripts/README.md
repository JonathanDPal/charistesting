I used these four scripts to help me out with keeping 
everything organized properly.

check.py: This script searched through a top level 
output directory indicated by command line argument, 
then indicated how many parameter combinations had
been completed and how many were left, based on hard
coded values for each parameter.

injection_locations.py: This script will take in some
hard coded suggestions for position angles and 
seperations for injected planets, along with an
observation datacube and the locations of its science
target(s). Then it will determine whether everything
is at least 2 FWHMs apart or not. If not, then 
locations of planet injections should be tweaked. 

nextparams.py: This script searched through a top level
directory, then suggested a set of parameter 
combinations which had not yet been completed. Since
I wasn't running all parameter combinations in one go
on the first couple sets of observations, this helped
me to finish filling out the parameter space.

standardize_filenames.py: This script was useful to 
me because slight shifts occurred in the filenames 
as parameter_test_infrastructure got some updates, 
but I needed to have everything standardized for 
analysis purposes. This script searches through a top 
level output directory and renames the files to be in 
standard format.