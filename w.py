import numpy as np

fake_fluxes = [5e-4, 5e-5, 5e-6, 1e-4, 1e-5, 1e-6]
fake_seps = [20, 40, 60]
fake_PAs=[19, 79, 139, 199, 259, 319]

a = [np.mean([fake_fluxes[0], fake_fluxes[3]]), np.mean([fake_fluxes[1], fake_fluxes[4]]), np.mean([fake_fluxes[2],
                                                                                                   fake_fluxes[5]])]

numgroups = len(fake_seps)
groupsize = int(len(fake_fluxes) / len(fake_seps))
groups = [[fake_fluxes[i * (groupsize + 1) + j] for i in range(groupsize)] for j in range(numgroups)]
fluxes = []
for i in range(numgroups):
    fxs = []
    for j in range(groupsize):
        fxs.append(groups[i][j])
    fluxes.append(np.mean(fxs))

print(fluxes)