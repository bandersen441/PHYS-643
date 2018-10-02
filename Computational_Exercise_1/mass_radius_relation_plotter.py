# This script plots the mass-radius relation for white dwarfs according to the data in the file
# pc_mass_radius.txt (which comes from ``hw1.cc'').
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc,rcParams
from matplotlib.ticker import ScalarFormatter
import matplotlib.patches as patches

rc('font',**{'family':'serif','serif':['Times']})
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16
rc('text', usetex=True)
fs=18
ms = 2
lw = 2

f = open('pc_mass_radius.txt')
lines = f.readlines()
f.close()

p_c = []
mass = []
radius = []
for line in lines[68:]:
    split = line.split('\t')
    p_c.append(float(split[0]))
    mass.append(float(split[1]))
    radius.append(float(split[2]))

# Plot the theoretical curve
M = np.linspace(0.0, 2.0, 5000)
R = 8.7 * 10**8 * M**(-1./3.) * (1 - (M / 1.43)**(4./3.))**0.5
plt.plot(M, R, color='green', linestyle='-', linewidth=2*lw, alpha = 0.5, zorder=0)
# Plot the R ~ M^(-1/3) low-mass approximation
R = 8.7 * 10**8 * M**(-1./3.)
plt.plot(M, R, color='blue', linestyle='-', linewidth=3*lw, alpha = 0.5, zorder=0)
# Plot the simulation data
plt.plot(mass, radius, color='black', linestyle='-', linewidth=lw, zorder=1)
# Plot a vertical line for the Chandrasehkar limit
plt.axvline(x=1.44, color='red', linestyle='--', linewidth=lw, zorder=0)
plt.gca().set_xlabel(r'Mass [M$_{\odot}$]',fontsize=fs,fontweight='bold')
plt.gca().set_ylabel(r'Radius [cm]',fontsize=fs,fontweight='bold')
plt.gca().set_yscale('log')
plt.gca().set_xlim([0.0,1.45])
plt.gca().set_ylim([2*10**7,6*10**9])
plt.tight_layout()
plt.savefig('mass_radius_relation.pdf', format='pdf', dpi=1000)
