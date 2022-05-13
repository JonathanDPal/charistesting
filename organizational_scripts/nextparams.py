import sys
from glob import glob
from valuefinder import valuefinder


def listset(alist):
    return list(set(alist))


direc = sys.argv[1]
klip_outputs = glob(f'{direc}/klipped_cubes_Wfakes/*.fits')
completed = [valuefinder(file, 'all') for file in klip_outputs]

fullset = []
for a in listset([pms[0] for pms in completed]):
    for s in listset([pms[1] for pms in completed]):
        for m in listset([pms[2] for pms in completed]):
            for p in listset([pms[3] for pms in completed]):
                for n in listset([pms[4] for pms in completed]):
                    for c in listset([pms[5] for pms in completed]):
                        for h in listset([pms[6] for pms in completed]):
                            if not [a, s, m, p, n, c, h] in completed:
                                fullset.append([a, s, m, p, n, c, h])

annuli = []
subsections = []
movement = []
spectrum = [None]
numbasis = []
smooth = []
highpass = []

first_time = True

for line in fullset:
    ann, sbs, mov, spec, nb, cs, hp = line
    if first_time:
        annuli.append(ann)
        subsections.append(sbs)
        movement.append(mov)
        numbasis.append(nb)
        smooth.append(cs)
        highpass.append(hp)
        first_time = False
    else:
        if ann not in annuli:
            good_to_go = True
            for S in subsections:
                for M in movement:
                    for N in numbasis:
                        for C in smooth:
                            for H in highpass:
                                if good_to_go and not [ann, S, M, None, N, C, H] in fullset:
                                    good_to_go = False
            if good_to_go:
                annuli.append(ann)
        if sbs not in subsections:
            good_to_go = True
            for A in annuli:
                for M in movement:
                    for N in numbasis:
                        for C in smooth:
                            for H in highpass:
                                if good_to_go and not [A, sbs, M, None, N, C, H] in fullset:
                                    good_to_go = False
            if good_to_go:
                subsections.append(sbs)
        if mov not in movement:
            good_to_go = True
            for A in annuli:
                for S in subsections:
                    for N in numbasis:
                        for C in smooth:
                            for H in highpass:
                                if good_to_go and not [A, S, mov, None, N, C, H] in fullset:
                                    good_to_go = False
            if good_to_go:
                movement.append(mov)
        if nb not in numbasis:
            good_to_go = True
            for A in annuli:
                for S in subsections:
                    for M in movement:
                        for C in smooth:
                            for H in highpass:
                                if good_to_go and not [A, S, M, None, nb, C, H] in fullset:
                                    good_to_go = False
            if good_to_go:
                numbasis.append(nb)
        if cs not in smooth:
            good_to_go = True
            for A in annuli:
                for S in subsections:
                    for M in movement:
                        for N in numbasis:
                            for H in highpass:
                                if good_to_go and not [A, S, M, None, N, cs, H] in fullset:
                                    good_to_go = False
            if good_to_go:
                smooth.append(cs)
        if hp not in highpass:
            good_to_go = True
            for A in annuli:
                for S in subsections:
                    for M in movement:
                        for N in numbasis:
                            for C in smooth:
                                if good_to_go and not [A, S, M, None, N, C, hp] in fullset:
                                    good_to_go = False
            if good_to_go:
                highpass.append(hp)

with open(f'{direc}_somenewparams.txt', 'w') as file:
    file.write(f'\nannuli = {annuli}')
    file.write(f'\nsubsections = {subsections}')
    file.write(f'\nmovement = {movement}')
    file.write(f'\nspectrum = {spectrum}')
    file.write(f'\nnumbasis = {numbasis}')
    file.write(f'\ncorr_smooth = {smooth}')
    file.write(f'\nhighpass = {highpass}')
