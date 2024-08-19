import numpy as np
from funcs import Discontinuous, Smooth
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
file="TestData/CjAAA10_o2.txt"
file="TestData/BSA_o2.txt"
file="TestData/blank_o2.txt"
data=np.loadtxt(file)
time=data[:,0]
current=data[:,1]
Smooth(time, current)

plt.show()
