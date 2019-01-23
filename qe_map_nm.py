# This code scan every grid moints in mesh
# This code also reomve the ions, once counted in a grid
# This code Do Not use adaptive meshing
# This code should show Normalized/raw QE on a grid


'''
Created on January 21, 2019
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

data= genfromtxt('laser_vs_anode_off/ionsmap_rc_6mm_laser_offset.txt',delimiter='',dtype=None, names=True)

e0_p=938.27e6

fig = plt.figure()

x=np.array(data['x'])
y=np.array(data['y'])
bz=np.array(data['Bz'])

Ez = e0_p*((1-bz**2)**(-0.5)-1)
Ezk = Ez/1000        # ion Energy in KeV

print(" Total Number of ion:", len(Ezk))

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
t = SY
print("t:", t)

#---- Creating Mesh Size-------
mesh = 0.00025

print("mesh size: ", mesh)

# -- Gird pints are created wth mesh value----
X = np.arange(-0.015, 0.015+mesh, mesh )
Y = np.arange(-0.015, 0.015+mesh, mesh)

#--- This XX, and YY is used later  for loop calculation
XX=X
YY=Y

#---- we used X,Y for grid, and XX, YY are left for loop calculation
X, Y = np.meshgrid(X, Y)
R= (X**2 + Y**2)*0

print("-- Code is Running; we will count till:--", len(XX)-1)


#-------------------------------------------------------------
# This function will check how many total ion points are
# nearby to a particular Grid Point. Once a ion is positioned at
# a certain grid, it wil not be used for later calcuation.
# it also removes ion, onc counted in acertain grid
def func_grid_check_tp(X_grid, Y_grid):
    points = 0
    remove_array =[]
    global x, y, t, mesh

    for k in range(len(x)):
        if abs(X_grid-x[k])<(mesh/2) and abs(Y_grid-y[k])<(mesh/2):
            points = points + 1

            remove_array.append(k)

    x = np.delete(x, remove_array)
    y = np.delete(y, remove_array)
    t = np.delete(t, remove_array)

    return points
#----------------------------------------------------------

#-----------------------------------------------------
# Thas part scan every grid points
for i in range(len(XX)):
    if len(t) == 0:
        break
    print("i : ", i)
    for j in range(len(YY)):
        R[j,i] = func_grid_check_tp(XX[i], YY[j])
#----------------------------------------------------

# -- Reverse Normalization ------------------------
# -- if omitted, will give raw QE value  ----------
Rmax = R.max()
Rmin = R.min()

print("-----------\nZmax:", R.max())
print("Zmin:", R.min())

for i in range(len(XX)):
    for j in range(len(YY)):
        R[j,i] = ((Rmax - R[j,i])/(Rmax - Rmin))
#-----------------------------------------------

Z=R
rcathode=0.013
plt.axes()

circle = plt.Circle((0, 0), radius=0.013, color='b',fill=False)
plt.gca().add_patch(circle)

plt.pcolor(X, Y, Z, cmap='gnuplot', clip_path=circle, clip_on=True)
#plt.pcolor(X, Y, Z, cmap='gnuplot2', clip_path=circle, clip_on=True)
plt.xlabel('X [m]')
plt.ylabel('Y [m]')

cb=plt.colorbar()
cb.set_label('QE/grid')

#plt.savefig('laser_vs_anode_off/ionsmap_rc_6mm_laser_offset_raw_qe_map_fine.png', dpi = 300)   # save the figure to file
#plt.close(fig)

# it prints total elapsed time in seconds
print("--- %s seconds ---" % (time.time() - start_time))

plt.show()


