# This is the plot for Sputterting Yield per 1000 ion
# This eq.s are obtained analytically by comparing with txt data set in gnuplot

from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
from numpy import genfromtxt
from traitlets.traitlets import add_article

fig = plt.figure()
ax1 = fig.add_subplot(111)



mesh = 0.1

X = np.arange(0, 350+mesh , mesh )
Y = np.arange(0, 350+mesh , mesh )
print("len(X):", len(X))



print("###########################")

for i in range(len(X)):
    if 0.2<X[i]<=1.1:
        Y[i]=36*(X[i]-0.2)
    elif 1.1<X[i]<=100:
        Y[i]=3.2+42*math.exp(-0.06*X[i])
    elif 100<X[i]<=350:
        Y[i]=-0.0044*X[i]+3.76
    else:
        Y[i]=0


#Y=Y/1000
plt.xlabel('Ion Energy[KeV]')
plt.ylabel('Cs Sputtering Yield [atom/ion] / 1000')
plt.title('Cs sputtering Yield with respect to H2+ ion Energy')


SRIM= genfromtxt('SRIM_data.txt',delimiter='',dtype=None, names=True)


xx=np.array(SRIM['HV'])
yy=np.array(SRIM['SY'])


ax1.plot(X, Y, '-.', label='fitted')
ax1.scatter(xx,yy, s=10, c='r', marker="o", label='SRIM data')
plt.legend(loc='upper right');



#plt.scatter(X, Y)
#plt.plot(X, Y,'-.')
#plt.plot(xx, yy,'.')
#plt.savefig('SPY/SP_Yield_vs_H2+_Ion_Energy.png', dpi = 300)
plt.show()
