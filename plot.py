import numpy as np 
import matplotlib.pyplot as plt 

dat = np.loadtxt('solution.dat')
x = dat[:,0]
y = dat[:,1]

plt.plot(x,y)
plt.show()