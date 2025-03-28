import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

sqrt_pi = 1.7724538509055160272981674833411451827975494561223871282138077898529112845

acc = 1e-4
eps = 1e-4
gauss = lambda x: np.exp(-x**2)
integral = quad(gauss, -np.Infinity, np.Infinity, epsabs=acc, epsrel=acc, full_output=1)
value = integral[0]
error = integral[1]
info = integral[2]
N = info['neval']

with open("out.txt", "a") as f:
    print(file=f)
    print("For the scipy.integrate.quad routine, same absolute and relative error: ", file=f)
    print("Value: ", value, file=f)
    print("Estimated error: ", error, file=f)
    print("Actual error: ", np.abs(value - sqrt_pi), file=f)
    print("Number of evaluations: ", N, file=f)
    print("So quad seems to underestimate the error, and so has an unnecesarryily large number of calls", file=f)