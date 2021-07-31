This is a set of scripts which will produce a set of 
plots based upon the planet detections. In order for
this to work, you'll need to have this folder located
one tier inside your top level output directory. For
example, if you named your top level output directory
'HD1160', then you'll want this folder to be 
'HD1160/analysis'. Once you set that up, you can run
'bash analysis.sh', and it will do everything for you
automatically. 

##### By everything, here is what it will produce:
1) 2D Heatmaps showing the average SNR of injected 
   planets as a function of annuli, subsections, 
   and movement.
2) 2D Heatmaps showing the average maximum SNR of 
   anything in the image as a function of annuli,
   subsections, and movement.
3) ROC curves showing the relationship between 
   true positive and false positive detections as 
   a function of number of KL modes, smoothing,
   and highpass.
   
Each one of these will be in a seperate subdirectory
within "{your object name}/analysis".
