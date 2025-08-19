import numpy as np
import matplotlib.pyplot as plt
nk=17
nw=100
data=np.loadtxt('Gpade.dat')
x=data[:,0].reshape(nk,nw)
y=data[:,1].reshape(nk,nw)
z=data[:,2].reshape(nk,nw)
plt.contourf(x,y,z,200)
plt.jet()
plt.colorbar()
plt.show()
