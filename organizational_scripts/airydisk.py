from scipy.special import j1, jn_zeros
import numpy as np


def airydisk(x, y, amplitude=1.0, x0= 0.0, y0 = 0.0, radius=1.0):
    rz = jn_zeros(1, 1)[0] / np.pi
    r = np.pi * np.sqrt((x - x0) ** 2 + (y - y0) ** 2) / (radius / rz)
    output = amplitude * ((2 * j1(r) / r) ** 2)
    return output
