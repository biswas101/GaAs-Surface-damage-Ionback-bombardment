#This code shows Normalized Power Density(Energy/grid) on Cathode position
# it is normalized f(n, E)/grid
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

data= genfromtxt('laser_vs_anode_off/ionsmap_rc_6mm_laser_offset.txt',delimiter='',dtype=None, names=True)

e0_p=938.27e6

fig = plt.figure()

x=np.array(data['x'])
y=np.array(data['y'])
bz=np.array(data['Bz'])

Ez = e0_p*((1-bz**2)**(-0.5)-1)
Ezk = Ez/1000        # ion Energy in KeV

print("len(Ezk):", len(Ezk))

SY = np.arange(0, len(Ezk) , 1 )   # SP. Yield, will be edited later


#---- this converstion is not necessary, in this script--------
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
print("-- Code is Running; we will count till:--", len(XX)-1)

#-----------------------------------------------------------------
# This function will check how many total ion points are
# nearby to a particular Grid Point and corresponding Energy.
# Once a ion energy is calculated at a certain grid, it wil not
# be used for later calcuation.
def func_grid_check_tp(X_grid, Y_grid):
    energy = 0
    remove_array =[]
    global x, y, t, mesh

    for k in range(len(x)):
        if abs(X_grid-x[k])<(mesh/2) and abs(Y_grid-y[k])<(mesh/2):
            energy = energy + t[k]  # adding total energy in a grid point

            remove_array.append(k)

    x = np.delete(x, remove_array)
    y = np.delete(y, remove_array)
    t = np.delete(t, remove_array)

    return energy
#----------------------------------------------------------------

#-----------------------------------------------------
# Thas part scan every grid points in mesh
for i in range(len(XX)):
    if len(t) == 0:
        break
    print("i : ", i)
    for j in range(len(YY)):
        R[j,i] = func_grid_check_tp(XX[i], YY[j])
#----------------------------------------------------


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
rlaser= 0.0043

fig,ax=plt.subplots()
plt.pcolor(X, Y, Z, cmap='gnuplot2_r')

plt.xlabel('X [m]')
plt.ylabel('Y [m]')
#plt.text(-0.014,0.0135,'xoffset 0mm - energy density')

cb=plt.colorbar()
cb.set_label('Energy/grid')
cir=plt.Circle((0,0),rcathode,color='k',fill=False)
#cir_ls=plt.Circle((0,0),rlaser,color='r',ls='dashed',fill=False)      # anode offset
cir_ls=plt.Circle((0.006,0),rlaser,color='r',ls='dashed',fill=False) # 6mm laser offset
ax.set_xlim((-0.015,0.015))
ax.set_ylim((-0.015,0.015))
ax.add_artist(cir)
ax.add_artist(cir_ls)


#plt.savefig('laser_vs_anode_off/ionsmap_rc_6mm_laser_offset_ED_map_nm_tight_3.png', dpi = 300 ,bbox_inches='tight')
#plt.close(fig)

# it prints total elapsed time in seconds
print("--- %s seconds ---" % (time.time() - start_time))

plt.show()





