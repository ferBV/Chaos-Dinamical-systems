import numpy as np
import matplotlib.pyplot as plt
import math as m
import scipy as sp
import seaborn as sns
from scipy.fft import fft, ifft

x = np.loadtxt("/home/fer/Documents/Caos/lorenz/seriex.txt")
y = np.loadtxt("/home/fer/Documents/Caos/lorenz/seriey.txt")
z = np.loadtxt("/home/fer/Documents/Caos/lorenz/seriez.txt")
t = np.loadtxt("/home/fer/Documents/Caos/lorenz/seriet.txt")

fx = fft(x)
fy = fft(y)
fz = fft(z)

plt.plot(t,x)
# plt.plot(fy)
# plt.plot(fz)

x2 = np.loadtxt("/home/fer/Documents/Caos/lorenz/seriex2.txt")
y2 = np.loadtxt("/home/fer/Documents/Caos/lorenz/seriey2.txt")
z2 = np.loadtxt("/home/fer/Documents/Caos/lorenz/seriez2.txt")

dx = np.absolute(np.add(x,-x2))
dy = np.absolute(np.add(y,-y2))
dz = np.absolute(np.add(z,-z2))
#plt.plot(dx)
# plt.plot(dy)
# plt.plot(dz)

