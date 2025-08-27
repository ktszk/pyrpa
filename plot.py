import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
nk=17
nw=100
data=np.loadtxt('Gpade.dat')
#nk=32
#nw=32
#data=np.loadtxt('Fk_tr.dat')
#data=np.loadtxt('gap_33.dat')

x=data[:,0].reshape(nk,nw)
y=data[:,1].reshape(nk,nw)
z=data[:,2].reshape(nk,nw)

min_col=z.min()
max_col=z.max()
if(min_col*max_col<0):
    crange=max_col-min_col
    zero_color=(0.-min_col)/crange
    mycolor=colors.LinearSegmentedColormap.from_list('custom_cmap',[(0.,'#0000FF'),(zero_color,'#FFFFFF'),(1.,'#FF0000')])
else:
    if min_col>0:
        mycolor=colors.LinearSegmentedColormap.from_list('custom_cmap',[(0.,'#FFFFFF'),(1.,'#FF0000')])
    elif max_col<0:
        mycolor=colors.LinearSegmentedColormap.from_list('custom_cmap',[(0.,'#FF0000'),(1.,'#FFFFFF')])
#plt.contourf(x,y,z,200,cmap=mycolor)
plt.contourf(x,y,z,200)
#plt.jet()
plt.hot()
plt.colorbar()
plt.show()
