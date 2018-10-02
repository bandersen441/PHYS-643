# This script creates a plot for each density profile dataset rk4_pc_####.txt and outputs them as
# density_prof_####.txt where #### is picked to order the files in terms of increasing mass.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc,rcParams
from matplotlib.ticker import ScalarFormatter
import matplotlib.patches as patches
import glob

rc('font',**{'family':'serif','serif':['Times']})
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16
rc('text', usetex=True)
fs=18
ms = 2
lw = 2

# Get filenames of all density profiles
files = sorted(glob.glob("./density_profiles/*.txt"))
n = len(files)
print files

max_r = []
min_r = []
max_p = []
min_p = []

progress = 0
for filename in files:
    f = open(filename)
    lines = f.readlines()
    f.close()

    radius = []
    density = []
    mass = []
    for line in lines[1:]:
        split = line.split('\t')

        radius.append(float(split[0]))
        density.append(float(split[1]))
        mass.append(float(split[2]))

    # For cosmetics
    density[-1] = 1

    # Actually plot the profile
    plt.plot(radius, density, color='black', linestyle='-', linewidth=lw, zorder=1)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.gca().set_xlim([10**5,6*10**9])
    plt.gca().set_ylim([10**(0),3*10**(13)])
    plt.gca().set_ylabel(r'Density [g cm$^{-3}$]', fontsize=fs, fontweight='bold')
    plt.gca().set_xlabel(r'Radius [cm]', fontsize=fs, fontweight='bold')
    x0, xmax = plt.xlim()
    y0, ymax = plt.ylim()
    plt.text(xmax-0.975*xmax,ymax-0.9*ymax, 'M = ' + "{0:.4f}".format(mass[-1]) + r' M$_{\odot}$', color='black',fontsize=fs)
    plt.subplots_adjust(left=0.15)
    plt.subplots_adjust(bottom=0.15)
    plt.savefig('./density_profiles/plots/log/density_profile_' + str(progress).zfill(4) + '.pdf', format='pdf', dpi = 300)
    plt.clf()
    progress = progress + 1
    print "Progress: " + str(float(progress)/n*100) + "%..."
