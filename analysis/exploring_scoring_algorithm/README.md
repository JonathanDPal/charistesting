Here I have three older scripts that I used for 
exploring the way that my scoring algorithm works. 
'compile_overall_scores.py' was a utility script which
took the overall scores which were generated as a weighted
average of the contrast score and SNR score for particular
parameters, and synthesized it across different star systems.
This overall score helped me to gain an initial feel for the
parameter space.
'scoreVtargetsnr.py' generates a CSV file with the contrast/snr
scores in one column and the sum of the science targets'
SNR values in the other. This allowed me to understand
how well the contrast & snr scores predict the actual
science targets' SNR. It takes one command line argument:
either 'snr' or 'contrast', depending on which scores you
are looking at. 'snr_contrast_scatter.py' will give you a 
plot with the SNR scores on the x-axis and the contrast 
scores on the y-axis. This helped me understand the interplay
between the two scores better.