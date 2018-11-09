# Bridget Andersen, 11/5/18
#    PHYS 643, Computational Exercise #3
#    This script contains all the code needed to derive the frequecies and eigenvalues
#    of the p-modes and g-modes of the Sun given a profile file from the 1M_pre_ms_to_wd
#    test suite in the MESA stellar evolution code.
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc,rcParams
from matplotlib.ticker import ScalarFormatter
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import sys

rc('text', usetex=True)
fs=18
generate_plots = False

# Potentially useful constants
R_sun = 6.96 * 10**(10)

# This function creates some basic plots of the initial MESA values
def plot_values(zone, R, p, P, g, H, cs, N2, identifier="init"):
    print "Plot values for " + identifier
    
    print "Plotting the density vs radius"
    plt.plot(R / R_sun, p, lw=1.5, color='r')
    plt.xlabel(r"r ($R_{\odot}$)")
    plt.ylabel(r"$\rho$ (g/cm$^3$)")
    plt.gca().set_yscale('log')
    plt.gca().set_xlim([0.,1.1])
    plt.savefig('p_R_' + identifier + '.pdf', format='pdf', dpi=1000, bbox_inches="tight", pad_inches=0.2)
    plt.clf()

    print "Plotting the pressure vs radius"
    plt.plot(R / R_sun, P * 10**(-9), lw=1.5, color='r')
    plt.xlabel(r"r ($R_{\odot}$)")
    plt.ylabel(r"P (GPa)")
    plt.gca().set_yscale('log')
    plt.gca().set_xlim([0.,1.1])
    plt.savefig('press_R_' + identifier + '.pdf', format='pdf', dpi=1000, bbox_inches="tight", pad_inches=0.2)
    plt.clf()

    print "Plotting the gravitational acceleration vs radius"
    plt.plot(R / R_sun, g, lw=1.5, color='r')
    plt.xlabel(r"r ($R_{\odot}$)")
    plt.ylabel(r"g (cm/s$^2$)")
    plt.gca().set_xlim([0.,1.1])
    plt.savefig('g_R_' + identifier + '.pdf', format='pdf', dpi=1000, bbox_inches="tight", pad_inches=0.2)
    plt.clf()

    print "Plotting the scale height vs radius"
    plt.plot(R / R_sun, H, lw=1.5, color='r')
    plt.xlabel(r"r ($R_{\odot}$)")
    plt.ylabel(r"H (cm)")
    plt.gca().set_yscale('log')
    plt.gca().set_xlim([0.,1.1])
    plt.savefig('H_R_' + identifier + '.pdf', format='pdf', dpi=1000, bbox_inches="tight", pad_inches=0.2)
    plt.clf()

    print "Plotting the sound speed vs radius"
    plt.plot(R / R_sun, cs, lw=1.5, color='r')
    plt.xlabel(r"r ($R_{\odot}$)")
    plt.ylabel(r"c$_s$ (cm/s)")
    plt.gca().set_xlim([0.,1.1])
    plt.savefig('cs_R_' + identifier + '.pdf', format='pdf', dpi=1000, bbox_inches="tight", pad_inches=0.2)
    plt.clf()

    print "Plotting the Brunt-Vaisala frequency vs radius"
    plt.plot(R / R_sun, N2, lw=1.5, color='r')
    plt.xlabel(r"r ($R_{\odot}$)")
    plt.ylabel(r"N$^2$")
    plt.gca().set_xlim([0.,1.1])
    plt.savefig('N2_R_' + identifier + '.pdf', format='pdf', dpi=1000, bbox_inches="tight", pad_inches=0.2)
    plt.clf()

# This function reads in the given MESA profile and separates it into
# arrays for each parameter that we will use in our p-modes and g-modes
def read_mesa_profile(filename='./solar_model.dat'):
    # Read in all lines
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    # Separate parameter header from the rest of the data and make the data a numpy
    # matrix that can be sliced
    params = lines[5].split()
    data = np.array([line.split() for line in lines[6:]])

    # Get the column index of each parameter
    idx_zone = params.index('zone')
    idx_logR = params.index('logR')
    idx_logp = params.index('logRho')
    idx_logP = params.index('logP')
    idx_logg = params.index('log_g')
    idx_H = params.index('pressure_scale_height')
    idx_cs = params.index('csound')
    idx_N2 = params.index('brunt_N2')
    # Make arrays for each parameter
    zone = np.array(data[:,idx_zone], dtype='int')
    R = np.power(10., np.array(data[:,idx_logR], dtype='double')) * R_sun
    p = np.power(10., np.array(data[:,idx_logp], dtype='double'))
    P = np.power(10., np.array(data[:,idx_logP], dtype='double'))
    g = np.power(10., np.array(data[:,idx_logg], dtype='double'))
    H = np.array(data[:,idx_H], dtype='double') * R_sun
    cs = np.array(data[:,idx_cs], dtype='double')
    N2 = np.array(data[:,idx_N2], dtype='double')

    # Change zone=1 from the surface to the center (makes more intuitive sense
    # once we get to integrating)
    R = R[::-1]
    p = p[::-1]
    P = P[::-1]
    g = g[::-1]
    H = H[::-1]
    cs = cs[::-1]
    N2 = N2[::-1]

    # Optionally plot parameters versus R
    if generate_plots:
        plot_values(zone, R, p, P, g, H, cs, N2)
    
    return zone, R, p, P, g, H, cs, N2

# Takes in all the relevant parameters, and outputs interpolated arrays
def interpolate_values(zone, R, p, P, g, H, cs, N2, n):
    zone = np.arange(1, n + 1)

    # Create interpolators
    interp_kind = 'linear'
    f_p = interp1d(R, p, kind=interp_kind)
    f_P = interp1d(R, P, kind=interp_kind)
    f_g = interp1d(R, g, kind=interp_kind)
    f_H = interp1d(R, H, kind=interp_kind)
    f_cs = interp1d(R, cs, kind=interp_kind)
    f_N2 = interp1d(R, N2, kind=interp_kind)
    # Create new interpolated values of R
    R = np.linspace(np.min(R), np.max(R), num=n)
    p = f_p(R)
    P = f_P(R)
    g = f_g(R)
    H = f_H(R)
    cs = f_cs(R)
    N2 = f_N2(R)

    # Optionally plot parameters versus R
    if generate_plots:
        plot_values(zone, R, p, P, g, H, cs, N2, identifier="interp")
    
    return zone, R, p, P, g, H, cs, N2

# Given an l and w (p0 and R0 too), determine the initial Er and dP
def initial_conditions(l, w, p0, R0):
    if l == 0:
        Er = 1.
        dP = 0.
    else:
        Er = 1.
        dP = p0 * w**2 * R0 * Er / l
    return Er, dP

# Compute the quantity d(Er)/dr
def dEr_dr(g, cs, R, Er, dP, p, l, w):
    return (g / cs**2 - 2. / R) * Er - (dP / p) * (1. / cs**2 - l * (l + 1.) / (w**2 * R**2))

# Compute the quantity d(dP)/dr
def ddP_dr(g, dP, cs, p, w, N2, Er):
    return -g * dP / cs**2 + p * (w**2 - N2) * Er

# Complete the Euler integration method
def euler_integration(R, p, P, g, H, cs, N2, l, w):
    # Get the initial values of Er and dP
    Er0, dP0 = initial_conditions(l, w, p[0], R[0])
    Er = [Er0]
    dP = [dP0]
    
    # Calculate the stepsize (will be determined by earlier interpolation)
    stepsize = R[1] - R[0]

    Er = list(np.copy(Er))
    dP = list(np.copy(dP))

    # Integrate from the center of the star and stop when we reach radius <= R_sun
    i = 0
    while R[i] <= R_sun:
        next_Er = Er[i] + dEr_dr(g[i], cs[i], R[i], Er[i], dP[i], p[i], l, w) * stepsize
        next_dP = dP[i] + ddP_dr(g[i], dP[i], cs[i], p[i], w, N2[i], Er[i]) * stepsize
        Er.append(next_Er)
        dP.append(next_dP)
        i = i + 1

    if generate_plots:
        ii = 0
        # Generate plots of eigenfunctions
        R_temp = R[:len(Er)]
        p_temp = p[:len(Er)]
        print "Plotting the Er vs radius for v=" + str(w * 10**(6) / (2. * np.pi))
        plt.gca().axhline(y=0, linewidth=1., color='black')
        plt.plot(R_temp / R_sun, R_temp * p_temp**(0.5) * Er, lw=1., color='r')
        plt.xlabel(r"r ($R_{\odot}$)")
        plt.ylabel(r"$r \rho^{1/2} \xi_r$")
        plt.gca().set_xlim([0.,1.])
        plt.title(r"Eigenfunction for $\ell = $" + str(l) + ", $v = $" + str(w * 10**(6) / (2. * np.pi)) + " $\mu$Hz")
        plt.savefig('./eigenfunction_v' + '{0:.2f}'.format(w * 10**(6) / (2. * np.pi)) + '.png', format='png', dpi=500, bbox_inches="tight", pad_inches=0.2)
        plt.clf()

        print "Plotting the dP vs radius"
        plt.gca().axhline(y=0, linewidth=1., color='black')
        plt.plot(R_temp / R_sun, dP, lw=1., color='r')
        plt.xlabel(r"r ($R_{\odot}$)")
        plt.ylabel(r"$\delta$P (Pa)")
        plt.title(r"$\ell = $" + str(l) + ", $v = $" + str(w * 10**(6) / (2. * np.pi)) + " $\mu$Hz")
        # plt.gca().set_yscale('log')
        plt.gca().set_xlim([0.,1.1])
        plt.savefig('dP_R_v' + str(ii).zfill(4) + '.png', format='png', dpi=1000, bbox_inches="tight", pad_inches=0.2)
        plt.clf()
    
    return Er, dP

if __name__ == '__main__':
    # Read in parameters that we will use in our Euler method
    zone, R, p, P, g, H, cs, N2 = read_mesa_profile()
    # Interpolate to n values
    n = 10 * len(zone) # 10 times the input length
    zone, R, p, P, g, H, cs, N2 = interpolate_values(zone, R, p, P, g, H, cs, N2, n)
    
    # Set the frequencies that we will iterate through
    # v = 100. * 10**(-6)
    vs = np.logspace(np.log10(40 * 10**(-6)), np.log10(4000 * 10**(-6)), 1000)
    ws = 2 * np.pi * vs
    # Set the l-values that we will iterate through
    # l = 1
    ls = [1]
    for l in ls:
        print "For l=" + str(l)
        # Define a short functions to use with brentq root finding
        def calc_bc(w):
            # Compute Euler integration and return the eigenvectors in Er and dP
            Er, dP = euler_integration(R, p, P, g, H, cs, N2, l, w)
            n = len(Er)
            # Return the boundary condition at r=R
            return dP[-1] / P[n-1] + Er[-1] / H[n-1]
        
        # Iterate through all the frequencies
        bc_arr = []
        for w in ws:
            # Find the value of the r=R boundary condition for the given frequency
            bc_arr.append(calc_bc(w))
        
        # Find all the roots for a given l
        # Find the indices when the values go from positive to negative or vice-versa
        sign_changes = np.where(np.diff(np.signbit(bc_arr)))[0]
        # Use scipy.brentq to find the root at each sign-changing interval
        roots = np.array([brentq(calc_bc, ws[index], ws[index + 1]) for index in sign_changes])
        n_roots = len(roots)
        print n_roots
        if generate_plots:
            print "Plotting the boundary condition vs w for l=" + str(l)
            plt.gca().axhline(y=0, linewidth=1., color='black')
            plt.plot(ws, bc_arr, lw=1., color='r')
            plt.xlabel(r"$\omega$")
            plt.ylabel(r"Boundary Condition")
            plt.scatter(roots, np.zeros(n_roots), c='blue', s=3, marker='o')
            plt.title(r"Boundary Condition and roots for $\ell = $" + str(l))
            plt.gca().set_xscale('log')
            plt.savefig('bc_w_l' + str(l) + '.pdf', format='pdf', dpi=1000, bbox_inches="tight", pad_inches=0.2)
            plt.clf()

        # For each eigenfunction, calculate the number of radial nodes
        n_nodes = []
        print roots
        colors = plt.cm.jet(np.linspace(0,1,n_roots))[::-1]
        plot = np.arange(0, n_roots, 5)
        for i in range(n_roots):
            Er, dP = euler_integration(R, p, P, g, H, cs, N2, l, roots[i])
            # Find the indices when the values go from positive to negative or vice-versa
            sign_changes = np.where(np.diff(np.signbit(Er)))[0]
            # Number of nodes is just the number of sign changes
            n_nodes.append(len(sign_changes))

            R_temp = R[:len(Er)]
            p_temp = p[:len(Er)]
            cs_temp = cs[:len(Er)]
            if i in plot:
                quant = R_temp * p_temp**(0.5) * Er
                plt.plot(R_temp / R_sun, quant / np.max(quant), lw=1., color=colors[i], alpha=0.5)
                plt.xlabel(r"r ($R_{\odot}$)")
                plt.ylabel(r"Normalized $r \rho^{1/2} \xi_r$")
        N2_temp = N2[:len(Er)]
        cs_temp = cs[:len(Er)]
        plt.plot(R_temp / R_sun, N2_temp / np.max(N2_temp), lw=1., linestyle='-.', color='black', alpha=1.)
        plt.gca().axhline(y=0, linewidth=1., color='black')
        plt.gca().set_xlim([0.,1.])
        plt.title(r"Multiple Eigenfunctions for $\ell = $" + str(l))
        plt.savefig('./mult_eigenfunctions.png', format='png', dpi=500, bbox_inches="tight", pad_inches=0.2)
        plt.clf()
        
        if generate_plots:
            print "Plotting the number of radial nodes vs w for l=" + str(l)
            plt.plot(np.array(roots) * 10**(6) / (2. * np.pi), n_nodes, lw=1., color='r')
            plt.xlabel(r"v ($\mu$Hz)")
            plt.ylabel(r"Number of Radial Nodes")
            plt.title(r"Number of Radial Nodes for g and p Mode Frequencies, $\ell = $" + str(l))
            plt.gca().set_xscale('log')
            plt.savefig('num_radial_nodes_l' + str(l) + '.pdf', format='pdf', dpi=1000, bbox_inches="tight", pad_inches=0.2)
            plt.clf()
        
    
