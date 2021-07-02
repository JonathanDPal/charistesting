import numpy as np

fake_fluxes = [1e-3, 1e-4, 1e-5, 5e-4, 5e-5, 5e-6] # List of Float(s)
fake_seps = [20, 40, 60] # List of Integer(s) and/or Float(s)
fake_PAs=[19, 79, 139, 199, 259, 319] # List of Integer(s) and/or Float(s)

# retrieved_fluxes = []
# for sep in fake_seps:
#     fake_planet_fluxes = []
#     for pa in fake_PAs:
#         fake_flux = retrieve_planet_flux(frame, dataset_center, output_wcs, sep, pa, searchrad=7)
#         fake_planet_fluxes.append(fake_flux)
#     retrieved_fluxes.append(np.mean(fake_planet_fluxes))

numgroups = len(self.fake_seps)
groupsize = int(len(self.fake_fluxes) / len(self.fake_seps))
fluxgroups = [[self.fake_fluxes[i*numgroups + j] for i in range(groupsize)] for j in range(numgroups)]
fluxes = [np.mean(fluxgroups[i]) for i in range(numgroups)]

print(fluxes)