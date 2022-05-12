This is a set of scripts which will produce a set of 
plots based upon the planet detections. In order for
this to work, you'll need to have this folder located
one tier inside your top level output directory. For
example, if you named your top level output directory
'HD1160', then you'll want this folder (currently named
'old_plotting_stuff') to be 'HD1160/analysis'. Once you
set that up, you can run 'bash analysis.sh {weight1} 
{weight2}', and it will do everything for you 
automatically. You'll need to specify the two weights, 
where {weight1} is the weight given to the SNR score in 
computing the overall score and {weight2} is the weight 
given to the contrast score. If you don't specify these,
then analysis.sh will still run everything except for 
the computing the overall score. Your weights can be 
on any scale you want; the only thing that matters is 
the ratio between them (e.g. there is no difference 
between specifying 0.25 & 0.75 and specifying 1 & 3).

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
4) CSV files which have the scores for each set of 
   parameters in terms of contrast, SNR, and overall.
   
Each one of these will be in a seperate subdirectory
within "{your object name}/analysis", assuming you 
set up the way I specified above.
