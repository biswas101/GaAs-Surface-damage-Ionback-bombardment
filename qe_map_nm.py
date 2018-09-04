#This code shows Normalized QE on Cathode position
# actually Normalized QE/grid

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



#data= genfromtxt('omer_offset_anode/ions_on_cathode_tab_0mm.txt',delimiter='',dtype=None, names=True)
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

print("SY:", SY)


#-- This part mimic SP. Yield with respect to ion Energy[KeV]

for i in range(len(Ezk)):   # SP due to H2+ ion
    if 0.2<Ezk[i]<=1.1:
        SY[i]=36*(Ezk[i]-0.2)
    elif 1.1<Ezk[i]<=100:
        SY[i]=3.2+42*math.exp(-0.06*Ezk[i])
    elif 100<Ezk[i]<=350:
        SY[i]= -0.0044*Ezk[i] +3.8
    else:
        SY[i]=0

#-----------------------------------------------------------

#---- this variable converstion is not necessary--------
x= xx
y= yy
t = SY
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

print("number of XX points:", len(XX))
print("number of YY points:", len(YY))

#---- we used X,Y for grid, and XX, YY are left for loop calculation
X, Y = np.meshgrid(X, Y)
R= (X**2 + Y**2)*0

print("##-- Code is Running, Please be Patient--##")

#--- this part will compare every point of grid, then it will
#-   find if any ion is in the grid or not. If ion is located
#-   in the grid, then SP. Yield corresponding to that ion-
#-   energy is added on the grid.
#----- For more than one ion in a grid, their SP. Yield
#-     is added simply

for i in range(len(XX)):
    for j in range(len(YY)):
        for k in range(len(Ezk)):
            if abs(XX[i]-x[k])<(mesh/2) and abs(YY[j]-y[k])<(mesh/2):
                R[j,i]= R[j,i] + t[k]
            else:
                R[j,i]= R[j,i]


# it prints total elapsed time in seconds
print("--- %s seconds ---" % (time.time() - start_time))

# -- Reverse Normalization ------------------------
Rmax = R.max()
Rmin = R.min()

print("Zmax:", R.max())
print("Zmin:", R.min())


for i in range(len(XX)):
    for j in range(len(YY)):
        R[j,i] = ((Rmax - R[j,i])/(Rmax - Rmin))


Z=R
#---------------------------------------------

rcathode=0.013

plt.axes()

circle = plt.Circle((0, 0), radius=0.013, color='b',fill=False)
plt.gca().add_patch(circle)



plt.pcolor(X, Y, Z, cmap='gnuplot', clip_path=circle, clip_on=True)
#plt.pcolor(X, Y, Z, cmap='gnuplot2', clip_path=circle, clip_on=True)
#plt.pcolor(X, Y, Z, cmap='gnuplot2_r', clip_path=circle, clip_on=True)   # -r will  reverse the colorbar
plt.xlabel('X [m]')
plt.ylabel('Y [m]')
#plt.text(-0.014,0.0135,'laser offset : 6mm')

cb=plt.colorbar()
cb.set_label('QE/grid')

plt.savefig('ionsmap_3r_vp_vzcut/ionsmap_rc_vp_20cm_qe_map_coarse.png', dpi = 300)   # save the figure to file
#plt.close(fig)

plt.show()


