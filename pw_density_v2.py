#This code shows Normalized Energy Density on Cathode position
# it is f(n, E)
# actually Normalized Energy/grid
# this will not clip the cathode, rather draw a circle over the cathode

'''
Created on August 6, 2018
@author:  Jyoti Biswas
Email:jbiswas@bnl.gov
'''

from matplotlib import cm
from matplotlib.colors import LogNorm
from matplotlib.mlab import bivariate_normal
from scipy.constants import *
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
from numpy import genfromtxt
import time
from traitlets.traitlets import add_article
import time
start_time = time.time()
from pylab import figure


data= genfromtxt('ionsmap_3r_vp_vzcut/ionsmap_rc_vp_20cm.txt',delimiter='',dtype=None, names=True)


e0_p=938.27e6

fig = plt.figure()

xx=np.array(data['x'])
yy=np.array(data['y'])
bz=np.array(data['Bz'])

Ez = e0_p*((1-bz**2)**(-0.5)-1)
Ezk = Ez/1000        # ion Energy in KeV

print("len(Ezk):", len(Ezk))

SY = np.arange(0, len(Ezk) , 1 )   # SP. Yield, will be edited later


#---- this converstion is not necessary, in this script--------
x= xx
y= yy
#t = SY
t=Ezk
print("t:", t)
#---- Creating Mesh Size-------
mesh = 0.0005
#mesh = 0.0002      # finer mesh takes longer time to execute

# -- Gird pints are created wth mesh value----
X = np.arange(-0.015, 0.015+mesh , mesh )
Y = np.arange(-0.015, 0.015+mesh, mesh)

#--- This XX, and YY is used later  for loop calculation
XX=X
YY=Y

print("number of XX Grid points:", len(XX))
print("number of YY Grid points:", len(YY))

#---- we used X,Y for grid, and XX, YY are left for loop calculation
X, Y = np.meshgrid(X, Y)
R= (X**2 + Y**2)*0
print("##-- Code is Running, Please be Patient--##")

#--- this part will compare every point of grid, then it will
#-   find if any ion is in the grid or not. If ion is located
#-   in the grid, then SP. Yield corresponding to that ion-
#-   energy is added on the grid points.
#----- For more than one ion in a grid, their Energy
#-     is added simply

for i in range(len(XX)):
    for j in range(len(YY)):
        for k in range(len(Ezk)):
            if abs(XX[i]-x[k])<(mesh/2) and abs(YY[j]-y[k])<(mesh/2):
                R[j,i]= R[j,i] + t[k]
            else:
                R[j,i]= R[j,i]


# total elapsed time in seconds
print("--- %s seconds ---" % (time.time() - start_time))



# -- This is Regular Normalization, range [0:15] --
Rmax = R.max()
Rmin = R.min()

print("Zmax:", R.max())
print("Zmin:", R.min())

for i in range(len(XX)):
    for j in range(len(YY)):
        R[j,i] = 15*((R[j,i] -Rmin)/(Rmax - Rmin))

Z=R
#---------------------------------------------

rcathode=0.013

fig,ax=plt.subplots()
plt.pcolor(X, Y, Z, cmap='gnuplot2_r')
#plt.pcolor(X, Y, Z, cmap='gnuplot2', clip_path=circle, clip_on=True)

plt.xlabel('X [m]')
plt.ylabel('Y [m]')
#plt.text(-0.014,0.0135,'xoffset 0mm - energy density')

cb=plt.colorbar()
#cb.set_label('Energy/grid')
cir=plt.Circle((0,0),rcathode,color='k',fill=False)
ax.set_xlim((-0.015,0.015))
ax.set_ylim((-0.015,0.015))
ax.add_artist(cir)


plt.savefig('ionsmap_3r_vp_vzcut/ionsmap_rc_vp_20cm_ED_map_tight_3.png', dpi = 300 ,bbox_inches='tight')
#plt.close(fig)

plt.show()


