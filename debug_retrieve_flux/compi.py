import numpy as np

with open('r2s.txt') as file:
    lines = [line for line in file]

fulldata = {str(sep): {str(pa): [] for pa in [19, 79, 139, 199, 259, 319]} for sep in [20, 40, 60]}


def valuefinder(s, p):
    paramlengths = {'sep': 3, 'pa': 2}
    if p == 'klip':
        assert s[0] == '['
        for j in range(len(s)):
            if s[j] == ']':
                finalindex = j
                break
        params = s[1:finalindex].split(',')
        modifiedparams = list()
        modifiedparams.append(int(params[0]))
        modifiedparams.append(int(params[1]))
        modifiedparams.append(float(params[2]))
        if params[3] == 'None':
            modifiedparams.append(None)
        else:
            raise ValueError("something went wrong -- spectrum isn't None")
        modifiedparams.append(int(params[4]))
        modifiedparams.append(float(params[5]))
        if params[6] == 'True':
            modifiedparams.append(True)
        elif params[6] == 'False':
            modifiedparams.append(False)
        else:
            modifiedparams.append(float(params[6]))
        return modifiedparams
    elif p in ['sep', 'pa']:
        pl = paramlengths[p]
        for i in range(len(s)):
            if s[i: i + pl] == p:
                subset = s[i+pl+1:]
                break
        if p == 'sep':
            val = subset.split('_')[0]
        else:
            val = subset.split('|')[0]
        return val
    else:
        return float(s.split('|')[-1])


for line in lines:
    sep, pa, r2 = valuefinder(line, 'sep'), valuefinder(line, 'pa'), valuefinder(line, 'r2')
    if r2 > 1:
        print(line)
    fulldata[sep][pa].append(r2)

means = {str(sep): {str(pa): np.mean(fulldata[str(sep)][str(pa)]) for pa in [19, 79, 139, 199, 259, 319]}
         for sep in [20, 40, 60]}
mins = {str(sep): {str(pa): np.min(fulldata[str(sep)][str(pa)]) for pa in [19, 79, 139, 199, 259, 319]}
         for sep in [20, 40, 60]}
maxs = {str(sep): {str(pa): np.max(fulldata[str(sep)][str(pa)]) for pa in [19, 79, 139, 199, 259, 319]}
         for sep in [20, 40, 60]}
stds = {str(sep): {str(pa): np.std(fulldata[str(sep)][str(pa)]) for pa in [19, 79, 139, 199, 259, 319]}
         for sep in [20, 40, 60]}

keypairs = [(sep, pa) for sep in [20, 40, 60] for pa in [19, 79, 139, 199, 259, 319]]

with open('avgr2s.csv', 'w') as file:
    file.write('sep, pa, avg, min, max, std')
for keypair in keypairs:
    sep, pa = keypair
    with open('avgr2s.csv', 'a') as file:
        file.write(f'\n{sep}, {pa}, {means[str(sep)][str(pa)]}, '
                   f'{mins[str(sep)][str(pa)]}, {maxs[str(sep)][str(pa)]}, {stds[str(sep)][str(pa)]}')
