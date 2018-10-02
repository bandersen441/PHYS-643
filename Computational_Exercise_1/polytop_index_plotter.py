# This script calculates the polytropic index gamma as the slope of ln(P) ~ ln(p) and
# plots it as a function of mass and density. Data is taken from "pc_mass_radius.txt"
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
ln_P = []
ln_p = []
for line in lines[68:]:
    split = line.split('\t')
    p_c.append(float(split[0]))
    mass.append(float(split[1]))
    ln_P.append(float(split[-2]))
    ln_p.append(float(split[-1]))

# Calculate gamma as the slope of ln(P) ~ ln(p)
gamma = []
pc_mid = []
mass_mid = []
for i in range(1,len(ln_P)):
    g = (ln_P[i] - ln_P[i-1])/(ln_p[i] - ln_p[i-1])
    gamma.append(g)
    pc_mid.append((p_c[i] + p_c[i-1]) / 2.)
    mass_mid.append((mass[i] + mass[i-1]) / 2.)

# Plot gamma vs mass
plt.plot(mass_mid, gamma, color='black', linestyle='-', linewidth=lw, zorder=0)
#plt.axvline(x=mass[792], color='red', linestyle=':', linewidth=lw, zorder=1)
plt.gca().set_xlabel(r'Mass [M$_{\odot}$]',fontsize=fs,fontweight='bold')
plt.gca().set_ylabel(r'$\gamma$',fontsize=fs,fontweight='bold')
plt.gca().set_xscale('log')
#plt.gca().set_yscale('log')
yax_ticks = np.around(np.arange(4./3.,5./3.,1./9.), decimals=2)
plt.gca().yaxis.set_ticks(yax_ticks)
plt.gca().set_xlim([1.5*10**(-3),2])
plt.gca().set_ylim([yax_ticks[0],yax_ticks[-1]])
plt.subplots_adjust(left=0.15)
plt.subplots_adjust(bottom=0.15)
plt.savefig('gamma_mass.pdf', format='pdf', dpi=1000)
plt.clf()

# Plot gamma vs density
plt.plot(pc_mid, gamma, color='black', linestyle='-', linewidth=lw, zorder=1)
#plt.axvline(x=10**6 / 0.5, color='red', linestyle=':', linewidth=lw, zorder=0)
plt.gca().set_xlabel(r'Density [g cm$^{-3}$]',fontsize=fs,fontweight='bold')
plt.gca().set_ylabel(r'$\gamma$',fontsize=fs,fontweight='bold')
plt.gca().set_xscale('log')
#plt.gca().set_yscale('log')
yax_ticks = np.around(np.arange(4./3.,5./3.,1./9.), decimals=2)
plt.gca().yaxis.set_ticks(yax_ticks)
plt.gca().set_xlim([10**(1),10**(13)])
plt.gca().set_ylim([yax_ticks[0],yax_ticks[-1]])
plt.subplots_adjust(left=0.15)
plt.subplots_adjust(bottom=0.15)
plt.savefig('gamma_density.pdf', format='pdf', dpi=1000)
