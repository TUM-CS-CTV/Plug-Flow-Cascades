import numpy as np
import pylab
tf = 100
tau = np.linspace(0,tf, num = 200) 

m1 = 0.05
m2 = 0.075
c = 2
tau2 = np.linspace(0,tf/c, num = 100) 
tau3 = np.linspace(tf/c,tf, num = 100) 

D1 = 10**(-5)
D2 = 10**(-7)
kA  = 2000 
kB  = 2000 
L   = 10**(-4) 
A =  10**(-13) 
N_alpha_A = 2*10**(13) / 2
N_alpha_B = 2*10**(13) / 2
N_beta = 2*10**(13)
EA_alpha_A = 10
EB_alpha_B = 10
EA_beta = EA_alpha_A/2
EB_beta = EB_alpha_B/2

m1_alpha_A = np.sqrt(kA*EA_alpha_A/D1)
m1_beta = np.sqrt(kA*EA_beta/D1)

m2_alpha_B = np.sqrt(kB*EB_alpha_B/D2)
m2_beta = np.sqrt(kB*EB_beta/D2)

p1_alpha_A = A * D1 * N_alpha_A * m1_alpha_A * np.tanh(m1_alpha_A*L)

p1_beta = A * D1 * N_beta * m1_beta * np.tanh(m1_beta*L)

p2_alpha_B = A * D2 * N_alpha_B * m2_alpha_B * np.tanh(m2_alpha_B*L)

p2_beta = A * D2 * N_beta * m2_beta * np.tanh(m2_beta*L)

p3_beta = A * D1 * N_beta * m1_beta**2 * m2_beta**2 /(m1_beta**2-m2_beta**2) * (np.tanh(m2_beta*L)/m2_beta - np.tanh(m1_beta*L)/m1_beta)

YI = 1 - np.exp(-p1_alpha_A*tau) - p1_alpha_A/(p2_alpha_B-p1_alpha_A) * (np.exp(-p1_alpha_A*tau) - np.exp(-p2_alpha_B*tau))

YII = (1 - np.exp(-p1_alpha_A*tau))*(1 - np.exp(-p2_alpha_B*tau))

YIII = 1-np.exp(-p1_beta*tau)-(p1_beta-p3_beta)/(p2_beta-p1_beta)*(np.exp(-p1_beta*tau)-np.exp(-p2_beta*tau))

import matplotlib.pyplot as plt
import matplotlib as mpl
fig = plt.figure(figsize=(9,4.3))
ax = fig.gca()
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['ytick.major.width'] = 3
ax.plot(tau,YI ,'-.', color = 'black', linewidth=2, label='$\it{Y}$$^\mathrm{I}$')
ax.plot(tau,YII, color = 'black', linewidth=2, label='$\it{Y}$$^\mathrm{II}$')
ax.plot(tau,YIII ,'--', color = 'black', linewidth = 2, label='$\it{Y}$$^\mathrm{III}$')
plt.xlabel('$\it{τ}$$_{\mathrm{f}}$ / (min)', fontsize=25)
plt.ylabel('$\it{Y}$$^\mathrm{m}$ / (mol/mol)', fontsize=25)
ax.grid(False)
ax.legend(fontsize=25)

plt.xticks(fontsize=25)
plt.xticks([0, 25, 50, 75, 100])
plt.yticks(fontsize=25)
ax.spines["top"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
plt.savefig('Case3.PNG', bbox_inches='tight')
plt.savefig('Case3.pdf', bbox_inches='tight')
#plt.show()