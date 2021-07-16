annuli = []
subsections = []
movement = []
spectrum = [None]
numbasis = []
smooth = []
highpass = []


def paramstring(ANN, SBS, MOV, NB, CS, HP):
    return f'[{ANN}, {SBS}, {MOV}, {None}, {NB}, {CS}, {HP}]\n'


lines = []
with open('suggestionsHR8799.txt', 'r') as file:
    for line in file:
        lines.append(line)

done = False
i = -1
first_time = True
while not done and i < (len(lines) - 1):
    i += 1

    line = lines[i]
    line = line.replace(' ', '')

    if line[0] != '[':
        continue
    else:
        ann, sbs, mov, spec, nb, cs, hp = line[1:-2].split(',')
        ann, sbs, mov, nb, cs = int(ann), int(sbs), float(mov), int(nb), float(cs)
        if hp == 'True':
            hp = True
        elif hp == 'False':
            hp = False
        else:
            hp = float(hp)
        if first_time:
            annuli.append(ann)
            subsections.append(sbs)
            movement.append(mov)
            numbasis.append(nb)
            smooth.append(cs)
            highpass.append(hp)
            first_time = False
        else:
            added_something = False
            if ann not in annuli:
                good_to_go = True
                for S in subsections:
                    for M in movement:
                        for N in numbasis:
                            for C in smooth:
                                for H in highpass:
                                    if good_to_go and not paramstring(ann, S, M, N, C, H) in lines:
                                        good_to_go = False
                if good_to_go:
                    annuli.append(ann)
                    added_something = True
            if sbs not in subsections:
                good_to_go = True
                for A in annuli:
                    for M in movement:
                        for N in numbasis:
                            for C in smooth:
                                for H in highpass:
                                    if good_to_go and not paramstring(A, sbs, M, N, C, H) in lines:
                                        good_to_go = False
                if good_to_go:
                    subsections.append(sbs)
                    added_something = True
            if mov not in movement:
                good_to_go = True
                for A in annuli:
                    for S in subsections:
                        for N in numbasis:
                            for C in smooth:
                                for H in highpass:
                                    if good_to_go and not paramstring(A, S, mov, N, C, H) in lines:
                                        good_to_go = False
                if good_to_go:
                    movement.append(mov)
                    added_something = True
            if nb not in numbasis:
                good_to_go = True
                for A in annuli:
                    for S in subsections:
                        for M in movement:
                            for C in smooth:
                                for H in highpass:
                                    if good_to_go and not paramstring(A, S, M, nb, C, H) in lines:
                                        good_to_go = False
                if good_to_go:
                    numbasis.append(nb)
                    added_something = True
            if cs not in smooth:
                good_to_go = True
                for A in annuli:
                    for S in subsections:
                        for M in movement:
                            for N in numbasis:
                                for H in highpass:
                                    if good_to_go and not paramstring(A, S, M, N, cs, H) in lines:
                                        good_to_go = False
                if good_to_go:
                    smooth.append(cs)
                    added_something = True
            if hp not in highpass:
                good_to_go = True
                for A in annuli:
                    for S in subsections:
                        for M in movement:
                            for N in numbasis:
                                for C in smooth:
                                    if good_to_go and not paramstring(A, S, M, N, C, hp) in lines:
                                        good_to_go = False
                if good_to_go:
                    highpass.append(hp)
                    added_something = True

print(f'Annuli: {annuli}')
print(f'Subsections: {subsections}')
print(f'Movement: {movement}')
print(f'Spectrum: {spectrum}')
print(f'Numbasis: {numbasis}')
print(f'Smooth: {smooth}')
print(f'Highpass: {highpass}')

with open('somenewparams.txt', 'w') as file:
    file.write(f'\nAnnuli: {annuli}')
    file.write(f'\nSubsections: {subsections}')
    file.write(f'\nMovement: {movement}')
    file.write(f'\nSpectrum: {spectrum}')
    file.write(f'\nNumbasis: {numbasis}')
    file.write(f'\nSmooth: {smooth}')
    file.write(f'\nHighpass: {highpass}')
