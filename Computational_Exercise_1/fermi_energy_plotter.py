# This script plots the central Fermi energy as a function of mass according to the
# data in the file pc_mass_radius.txt (which comes from ``hw1.cc'')
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

m_e = 9.1093897 * 10.**(-28);
c = 2.99792458 * 10.**(10);

p_c = []
mass = []
radius = []
E_F = []
for line in lines[68:]:
    split = line.split('\t')
    p_c.append(float(split[0]))
    mass.append(float(split[1]))
    radius.append(float(split[2]))
    E_F.append(float(split[3]))
E_F = np.array(E_F)
# Plot the Fermi energy ratio vs mass
plt.plot(mass, E_F / (m_e * c**2), color='black', linestyle='-', linewidth=lw, zorder=0)
# Plot a vertical line for the degeneracy limit (mass when p_c ~ 10^6 / 0.5 g cm^-3)
plt.axvline(x=mass[792], color='red', linestyle=':', linewidth=lw, zorder=1)
print mass[792]
print p_c[792]
plt.gca().set_xlabel(r'Mass [M$_{\odot}$]',fontsize=fs,fontweight='bold')
plt.gca().set_ylabel(r'E$_F / m_e c^2$',fontsize=fs,fontweight='bold')
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().set_xlim([1.5*10**(-3),2])
plt.gca().set_ylim([10**(-4),1.5*10**2])
plt.subplots_adjust(left=0.15)
plt.subplots_adjust(bottom=0.15)
plt.savefig('fermi_energy_mass.pdf', format='pdf', dpi=1000)
plt.clf()

# Plot the Fermi energy ratio vs density
plt.plot(p_c, E_F / (m_e * c**2), color='black', linestyle='-', linewidth=lw, zorder=1)
plt.axvline(x=10**6 / 0.5, color='red', linestyle=':', linewidth=lw, zorder=0)
plt.gca().set_xlabel(r'Density [g cm$^{3-3}$]',fontsize=fs,fontweight='bold')
plt.gca().set_ylabel(r'E$_F / m_e c^2$',fontsize=fs,fontweight='bold')
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().set_xlim([10**(1),10**(13)])
plt.gca().set_ylim([10**(-4),1.5*10**2])
plt.subplots_adjust(left=0.15)
plt.subplots_adjust(bottom=0.15)
plt.savefig('fermi_energy_density.pdf', format='pdf', dpi=1000)
