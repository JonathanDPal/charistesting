I used these scripts to help me out with keeping 
everything organized properly.

Short Version:
-
"valuefinder" is incredibly useful for basically any sort
of analysis that you may perform. Aside from that, it's
unlikely that any of these scripts would be useful to you.

Full Version:
-

check.py: This script searched through a top level 
output directory indicated by command line argument, 
then indicated how many parameter combinations had
been completed with respect to KLIP and how many were
left, based on hard coded values for each parameter. 
It also creates a text file containing all the remaining
parameter combinations, where each parameter combination
is on one line, with individual parameters seperated by
commas.

check_contrast_detections.py: Same thing as check.py, but
with contrast and/or detections instead of KLIP runs. 

get_dataset_info.py: This script will generate a text file with
information about the CHARISData class that would be built out
of the FITS files associated with the observation on some star
system, which you specify as a command line argument. If you
aren't going to be running KLIP, then you can use this to build
a TestDataset class which has the necessary info to do
contrast/detections, without having to build an entire CHARISData
class like you normally would.

move_and_delete.py: I used this script when using
the ND CRC's clusters so that I could automate
the file transfer.

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

valuefinder.py: Has two functions which are used in some of
the above scripts. They are two of the most useful things in
this entire repository. Will be helpful to you if you are 
trying to run your own analyses which aren't covered by what
I have written.
