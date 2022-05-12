I made these to help with finding the best places
to inject fake planets.

check_injection_locations.py: This script will 
take in some hard coded suggestions for position angles 
and seperations for injected planets, along with an
observation datacube and the locations of its science
target(s). Then it will determine whether everything
is at least 2 FWHMs apart or not. If not, then 
locations of planet injections should be tweaked. 

optimal_injection_locations.py: This script is similar 
to the previous one, except that instead of checking whether 
your suggestions work, you will give it your suggestions
and some bounds, and then it will optimize the locations 
of the fake planets.
