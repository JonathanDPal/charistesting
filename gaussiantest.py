import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

np.random.seed(314)
X = np.linspace(0, 4, 20)
Y = np.linspace(0, 4, 20)
xvals = np.array(list(X) * 20)
yvals = np.array([[yval] * 20 for yval in Y]).flatten()

fwhm = 3.5
peak = 75


def gaussian(xy, peak, fwhm, offset):
    x, y = xy
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return peak * np.exp(-(x ** 2 + y ** 2) / (2 * sigma ** 2)) + offset


df = pd.DataFrame()
df['x'] = xvals
df['y'] = yvals
modgauss = [gaussian((x, y), peak, fwhm, 0) if ((x ** 2 + y ** 2 > 1) and (x ** 2 + y ** 2 < 4))
            else (20 + ((-1 ** np.random.randint(2)) * np.random.randint(5))) for x, y in zip(xvals, yvals)]
df['mod gaussian'] = modgauss
print([np.sqrt(x, y) for x, y in zip(df['x'], df['y'])])
raise KeyboardInterrupt
df['distance'] = [np.sqrt(x, y) for x, y in zip(df['x'], df['y'])]
df = df[df['distance'] < 4]

coordinates = (df['x'], df['y'])

optparams, covmatrix = curve_fit(gaussian, coordinates, df['mod gaussian'])

xpred = np.linspace(0, 4, 3000)
ypred = np.linspace(0, 4, 3000)
xpredvals = np.array(list(xpred) * 3000)
ypredvals = np.array([[yval] * 3000 for yval in ypred]).flatten()
predictions = [gaussian((x, y), optparams[0], optparams[1], optparams[2]) for x, y in zip(xpredvals, ypredvals)]
preddistances = [np.sqrt(x, y) for x, y in zip(xpredvals, ypredvals)]

plt.scatter(df['distance'], df['mod gaussian'])
plt.plot(preddistances, predictions)
plt.show()