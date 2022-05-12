Here I have two older scripts that I used for 
exploring the way that my scoring algorithm works. 
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