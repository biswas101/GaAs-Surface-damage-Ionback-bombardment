# This code make a rough scan using big mesh, and later make a deep scan on big meshes where ions are located.
# While deep sacn, This code reomve the ions, once counted in a small grid
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


#data= genfromtxt('omer_offset_anode/ions_on_cathode_tab_0mm.txt',delimiter='',dtype=None, names=True)
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
big_mesh = 0.001
mesh_factor = 4
mesh = big_mesh/mesh_factor  # ie. 0.000250 or, 0.0000625

print("small mesh size: ", mesh)

# -- Gird pints are created wth mesh value----
X = np.arange(-0.015, 0.015+mesh, mesh )
Y = np.arange(-0.015, 0.015+mesh, mesh)

#print("X grid points:", X)
#print("Y grid points:", Y)

#--- This XX, and YY is used later  for loop calculation
XX=X
YY=Y

#---- we used X,Y for grid, and XX, YY are left for loop calculation
X, Y = np.meshgrid(X, Y)
R= (X**2 + Y**2)*0

print("##-- Code is Running, Please be Patient--##")

#X_big = np.arange(-0.015, 0.015+big_mesh , big_mesh )
#Y_big = np.arange(-0.015, 0.015+big_mesh, big_mesh)

# big grid points are positioned in a way to match the small grid points scan area.
X_big = np.arange(-0.015+(mesh_factor-1)*(mesh/2), 0.015+big_mesh-(mesh_factor-1)*(mesh/2) , big_mesh )
Y_big = np.arange(-0.015+(mesh_factor-1)*(mesh/2), 0.015+big_mesh-(mesh_factor-1)*(mesh/2), big_mesh)


XX_b=X_big
YY_b=Y_big

X_big, Y_big = np.meshgrid(X_big, Y_big)
R_b= (X_big**2 + Y_big**2)*0

xgrid_array= []
ygrid_array= []


# This part make a rough scan and find in which grid contain ions
for n in range(len(XX_b)):  # n is the number of big X grid point, contain one or more ions
    print( n, end=":", flush=True)
    for l in range(len(YY_b)):
        for m in range(len(t)):   # len(Ezk)
            if abs(XX_b[n]-x[m])<(big_mesh/2) and abs(YY_b[l]-y[m])<(big_mesh/2):
                #R_b[l,n]= R_b[l,n] + t[m]
                R_b[l, n] = 1
                xgrid_array.append(n)
                ygrid_array.append(l)
                break


#print("xgrid_array:", xgrid_array)
#print("ygrid_array:", ygrid_array)


print("\nRaw scan completed. \nWe are doing deep scan now!!")
print("Len xgrid/ygrid_array:", len(xgrid_array))
#print("Len ygrid_array:", len(ygrid_array))


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


#----------------------------------------------------
# This part does, deep scan in big grid points, where ions are located
for nn in range(len(xgrid_array)):
    if len(t) == 0:
        break
    print("nn : ", nn)
    #print("xgrid_array[nn]:", xgrid_array[nn])
    #print("ygrid_array[nn]:", ygrid_array[nn])
    for ii in range(mesh_factor):
        for jj in range(mesh_factor):

            i = xgrid_array[nn]*mesh_factor + ii
            j = ygrid_array[nn]*mesh_factor + jj
            #print("i, j: ", i, j)
            R[j, i] = func_grid_check_tp(XX[i], YY[j])
#-----------------------------------------------------

#-----------------------------------------------------
# Thas part scan every grid points
#for i in range(len(XX)):
#    if len(t) == 0:
#        break
#    print("i : ", i)
#    for j in range(len(YY)):
#        R[j,i] = func_grid_check_tp(XX[i], YY[j])
#----------------------------------------------------


# -- Reverse Normalization ------------------------
# -- if omitted, will give raw QE value  ----------
#'''
Rmax = R.max()
Rmin = R.min()

print("-----------\nZmax:", R.max())
print("Zmin:", R.min())


for i in range(len(XX)):
    for j in range(len(YY)):
        R[j,i] = ((Rmax - R[j,i])/(Rmax - Rmin))
#-----------------------------------------------
#'''


Z=R
Z_b=R_b # in case, you may want to see the mig scan area

rcathode=0.013

plt.axes()

circle = plt.Circle((0, 0), radius=0.013, color='b',fill=False)
plt.gca().add_patch(circle)


#plt.pcolor(X_big, Y_big, Z_b, cmap='gnuplot_r', clip_on=True) # you can see the big scan area using this
plt.pcolor(X, Y, Z, cmap='gnuplot', clip_on=True)
#plt.pcolor(X, Y, Z, cmap='gnuplot_r', clip_path=circle, clip_on=True)
#plt.pcolor(X, Y, Z, cmap='gnuplot2', clip_path=circle, clip_on=True)
#plt.pcolor(X, Y, Z, cmap='gnuplot2_r', clip_path=circle, clip_on=True)   # -r will  reverse the colorbar
plt.xlabel('X [m]')
plt.ylabel('Y [m]')
#plt.text(-0.014,0.0135,'laser offset : 6mm')

cb=plt.colorbar()
cb.set_label('QE/grid')

#plt.savefig('laser_vs_anode_off/ionsmap_rc_6mm_laser_offset_raw_qe_map_fine.png', dpi = 300)   # save the figure to file
#plt.close(fig)

# it prints total elapsed time in seconds
print("--- %s seconds ---" % (time.time() - start_time))

plt.show()


