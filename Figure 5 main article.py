import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib as klp
 
xlist = np.logspace(-2.10,2.10, num = 500) 
ylist = np.logspace(-2.05,2.05, num = 500) 

l1, l2 = np.meshgrid(xlist, ylist) 


Y1z = 1 - np.exp(-l1) - l1/(l2-l1) * (np.exp(-l1) - np.exp(-l2))
Y2z = (1 - np.exp(-l1))*(1 - np.exp(-l2))
R = Y1z/Y2z


from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

cmap = plt.cm.Blues(np.linspace(0,1.1,35))
cmap = klp.colors.ListedColormap(cmap[8:,:-1])
plt.rcParams.update({'font.size': 20})
plt.rcParams['axes.linewidth'] = 1.4 

plt.rcParams["figure.figsize"] = [6.1, 5]
plt.rcParams['figure.dpi'] = 1000

levels = [0.2, 0.4, 0.6, 0.8]

fig,ax=plt.subplots(1,1)
cp = ax.contourf(l1, l2, R, [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], vmin = 0.0, vmax=1.0, cmap=cmap.reversed())
cp2 = ax.contour(l1, l2, Y2z,  levels=levels, colors = 'white', linestyles = 'dashed', linewidths = 1.5)
manual_locations = [(2*10**(1),1.5*10**(-1)),(2*10**(1),3*10**(-1)),(2*10**(1),1*10**(0)),(2*10**(1),4*10**(0))]
#manual_locations = [(5*10**(0),1.5*10**(-1)),(3*10**(1),4*10**(-1)),(5*10**(0),6*10**(-1)),(3*10**(1),4*10**(0))]


fig.colorbar(cp, label='$R^{\mathrm{I/II}}$ / (mol/mol)') # Add a colorbar to a plot

plt.yscale('log')
plt.xscale('log')


ax.yaxis.set_tick_params(width=1.4, length = 5)
ax.xaxis.set_tick_params(width=1.4, length = 5)

ax.tick_params(axis='y', which='minor', colors='black', left=True, length = 3, width=1.2)
ax.tick_params(axis='x', which='minor', colors='black', bottom=True, length = 3, width=1.2)
ax.tick_params(axis="x", direction="out")


locmin = klp.ticker.LogLocator(base=10.0,subs=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=8)
ax.xaxis.set_minor_locator(locmin)
ax.xaxis.set_minor_formatter(klp.ticker.NullFormatter())


ax.set_xlim(10**(-2), 10**(2))
ax.set_ylim(10**(-2), 10**(2))

plt.yticks([10**(-2), 10**(-1), 10**(0), 10**(1), 10**(2)])
plt.xticks([10**(-2), 10**(-1), 10**(0), 10**(1), 10**(2)])

ax.set_xlabel('$\mu_1$ / (-)')
ax.set_ylabel('$\mu_2$ / (-)')
ax.clabel(cp2, levels, inline=1, fmt='%1.1f', fontsize=20, manual=manual_locations)

plt.savefig('contour12.png', bbox_inches='tight')  
plt.savefig('contour12.pdf', bbox_inches='tight')   
plt.show()