import numpy as np
from funcs import ChronoGui
data=np.loadtxt("CjAAA10_o2.txt")
time=data[:,0]
current=data[:,1]
ChronoGui(time, current, dumpfile="Example.csv")