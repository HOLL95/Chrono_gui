import numpy as np
from funcs import Discontinuous, Smooth
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
file="TestData/chrono_Cf2.txt"
data=np.loadtxt(file)
time=data[:,0]
current=data[:,1]
Smooth(time, current)

plt.show()
