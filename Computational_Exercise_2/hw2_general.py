# Bridget Andersen, 10/20/18
#    PHYS 643, Computational Exercise #2
#    This script contains all the code for a 1D hydrodynamics simulation implemented
#    using a simple donor cell advection scheme. We consider the conservation of mass
#    and momentum.
import numpy as np
import math
import matplotlib.pyplot as plt
import argparse

# This function takes in the current arrays for f1 and v and returns the next f1 profile
def update_f1(f1, v):
    # Calculate the velocities at the lower boundary of each cell
    v_low = 0.5 * (v + np.roll(v, 1))
    # We create masks that will tell us how to calculate the flux at the lower boundary of each cell,
    # depending on the velocity at the lower boundary
    pos_mask = v_low > 0
    neg_mask = v_low <= 0
    # Determine the Js at the lower boundary of each cell according to masks
    J_low = np.zeros(len(f1))
    J_low[pos_mask] = (v_low * np.roll(f1, 1))[pos_mask]
    J_low[neg_mask] = (v_low * f1)[neg_mask]

    # Since the lower boundary flux of one cell is equal to the higher boundary flux of another
    # cell, we can simply write:
    J_high = np.roll(J_low, -1)

    # Finally calculate the new f1 values
    f1_new = f1 - (dt / dx) * (J_high - J_low)
    return f1_new

# This function takes in arrays for the future f1 values (previously calculated when f1 was updated, see update_f1()),
# the f1 values of the current step, the f2 values of the current step, and the velocity values of the current step
# and returns the next f2 profile
def update_f2(f1_new, f2, v):
    # The first part of the update is just the standard donor-cell advection as defined in update_f1()
    f2_new = update_f1(f2, v)
    # Now add in the pressure gradient term, using the second order approximation
    f2_new = f2_new - (dt / (2. * dx)) * cs**2 * (np.roll(f1_new, -1) - np.roll(f1_new, 1))
    return f2_new

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--num_x", nargs='?', type=int, default=100+2, help="The number of grid points (int).", required=False)
    parser.add_argument("--n_iter", nargs='?', type=int, default=1000, help="The number of iterations to complete (int).", required=False)
    parser.add_argument("--dt", nargs='?', type=float, default=0.01, help="The timestep (float).", required=False)
    parser.add_argument("--p_amp", nargs='?', type=float, default=1., help="The initial density amplitude (float).", required=False)
    parser.add_argument("--v_amp", nargs='?', type=float, default=1., help="The initial velocity amplitude (float).", required=False)
    parser.add_argument("--cs", nargs='?', type=float, default=1., help="The initial sound speed (float).", required=False)
    parser.add_argument("--pause", nargs='?', type=float, default=0.001, help="How long to pause between each frame (float).", required=False)
    parser.add_argument("--v_type", nargs='?', type=str, default='zeros', help="The profile of the initial velocity (str: constant, gaussian, sinusoid, zeros).", required=False)
    parser.add_argument("--p_type", nargs='?', type=str, default='gaussian', help="The profile of the initial density (str: gaussian, sinusoid).", required=False)
    parser.add_argument("--example_scaling", nargs='?', type=int, default=1, help="1 indicates to use the specific plot scaling chosen for the examples in this exercise write-up. 0 will allow automatic scaling (int: 1 or 0).", required=False)
    args = parser.parse_args()

    # Parse all the arguments
    num_x = args.num_x + 2
    n_iter = args.n_iter
    dt = args.dt
    p_amp = args.p_amp
    v_amp = args.v_amp
    cs = args.cs
    pause = args.pause
    v_type = args.v_type
    p_type = args.p_type
    example_scaling = args.example_scaling

    # Some initial values
    dx = 100. / float(num_x) # The size of each x grid point
    PI = math.pi
    w = 2 * PI / (num_x * dx / 4.) # The angular frequency of the initial waves
    
    # Set the initial conditions
    x = np.arange(num_x) * dx
    # For the velocity
    if v_type == 'constant':
        v = v_amp * np.ones(len(x))
    elif v_type == 'gaussian':
        v = v_amp * np.exp(-(x-50)**2/5.**2)
    elif v_type == 'sinusoid':
        v = v_amp * np.sin(0.5 * w * x)**2 + 0.01
        v_low = -abs(1.5*np.min(v))
        v_high = abs(1.5*np.max(v))
    elif v_type == 'zeros':
        v = np.zeros(len(x))
        v_low = -0.05
        v_high = 0.05
    # For the density
    if p_type == 'gaussian':
        f1 = p_amp * np.exp(-(x-50)**2/5.**2) + 1.
        f1_low = 1.
        f1_high = np.max(f1) + 0.01
    elif v_type == 'sinusoid':
        f1 = p_amp * np.sin(0.5 * w * x)**2 + 1
        f1_low = 0.
        f1_high = 3*np.max(f1)
    f2 = f1 * v
    
    # Set up the plots and plot the initial conditions
    plt.ion()
    fig = plt.figure(figsize=(8, 5))
    ax1 = plt.subplot(121)
    dat1, = ax1.plot(x[1:-1], f1[1:-1], marker='o', c='blue', markersize=2)
    if example_scaling == 1:
        ax1.set_xlim([0., 100.])
        ax1.set_ylim([f1_low, f1_high])
    ax1.set_xlabel('x grid')
    ax1.set_ylabel('Density')

    ax2 = plt.subplot(122)
    dat2, = ax2.plot(x[1:-1], v[1:-1], marker='o', c='red', markersize=2)
    if example_scaling == 1:
        ax2.set_xlim([0., 100.])
        ax2.set_ylim([v_low, v_high])
    ax2.set_xlabel('x grid')
    ax2.set_ylabel('Velocity')

    plt.tight_layout()
    plt.show()
    plt.pause(pause)
    
    # Iterate through each timestep
    for t in range(n_iter):
        # Update f1=mass density
        f1_new = update_f1(f1, v)

        # Update f2=momentum density
        f2_new = update_f2(f1_new, f2, v)
        # Get the new velocity vector
        v_new = f2 / f1

        # Note that periodic boundary conditions are already implemented because of our
        # use of numpy.roll()

        # plot new values and pause for viewing
        dat1.set_ydata(f1[1:-1])
        dat2.set_ydata(v[1:-1])
        plt.draw()
        plt.pause(pause)

        # Set the current values to the new values
        f1 = f1_new
        f2 = f2_new
        v = v_new
