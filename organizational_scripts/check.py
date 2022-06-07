import os
import sys
import numpy as np
from glob import glob
from valuefinder import valuefinder, paramvaluesfinder

direc = sys.argv[1]
klip_outputs = glob(f'{direc}/klipped_cubes_Wfakes/*.fits')
try:
    completed = [valuefinder(file, 'all') for file in klip_outputs]
except ValueError:  # if there's new stuff and we haven't gotten rid of collapsed KLIP files yet
    os.system(f'mv {direc}/klipped_cubes_Wfakes/*modes-all* {direc}/kl-all')
    completed = [valuefinder(file, 'all') for file in klip_outputs]

olddirec = os.getcwd()
os.chdir(direc)
try:
    annuli = paramvaluesfinder('annuli')
except UnboundLocalError:  # if the log file got replaced by some new janky thing
    os.chdir(olddirec)
    os.system(f'cp useful_logfile.txt {direc}/log.txt')
    os.chdir(direc)
    annuli = paramvaluesfinder('annuli')

subsections = paramvaluesfinder('subsections')
movement = paramvaluesfinder('movement')
numbasis = paramvaluesfinder('numbasis')
corr_smooth = paramvaluesfinder('corr_smooth')
highpass = paramvaluesfinder('highpass')

ce = list()
nc = list()
for ann in annuli:
    for sbs in subsections:
        for mov in movement:
            for spec in [None]:
                for nb in numbasis:
                    for cs in corr_smooth:
                        for hp in highpass:
                            if [ann, sbs, mov, spec, nb, cs, hp] in completed:
                                ce.append([ann, sbs, mov, spec, nb, cs, hp])
                            else:
                                nc.append([ann, sbs, mov, spec, nb, cs, hp])
print(len(ce))  # number completed
print(len(ne))  # number left

with open(f'remainingparams.txt', 'w') as file:
    for incompleteparams in nc:
        ann, sbs, mov, spec, nb, cs, hp = incompleteparams
        file.write(f'{ann},{sbs},{mov},{spec},{nb},{cs},{hp}\n')
if len(sys.argv) > 2 and sys.argv[2] == 'complete':
    with open(f'completedparams.txt', 'w') as file:
        for completeparams in ce:
            ann, sbs, mov, spec, nb, cs, hp = completeparams
            file.write(f'{ann},{sbs},{mov},{spec},{nb},{cs},{hp}\n')
