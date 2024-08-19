import numpy as np
from funcs import Discontinuous
data=np.loadtxt("TestData/blank chrono_2.txt")
time=data[:,0]
current=data[:,1]
Discontinuous(time, current, dumpfile="Example.csv")