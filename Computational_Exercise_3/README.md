# PHYS 643 - Computational Exercise 3

The goal of this computational exercise was to derive the frequecies and eigenvalues of the p-modes and g-modes of the Sun given a profile file from the 1M_pre_ms_to_wd test suite in the MESA stellar evolution code.

## Directory contents
* As per in-class request, I have not completed a write-up for this assignment
* **hw3.py** : the python file containing all of my code for this assignment. I use Euler integration and scipy.optimize.brentq() to find the eigenfunctions for l=1.
* **bc_w_l1.pdf** : A plot showing the boundary condition evaluated at r=R_sun for many frequencies. The blue dots along y=0 show the roots solved for using scipy.optimize.brentq()
* **mult_eigenfunctions.png** : A plot showing many of the eigenfunctions for l=1 in terms of their radial displacement. The g-modes are shown in a more reddish color (lower frequencies) while the p-modes are shown in a more blueish color (higher frequencies). The black dotted line shows the Brunt-Vaisala frequency squared.
* **num_radial_nodes_l1.pdf** : A plot showing the number of radial nodes compared to the frequencies of the g and p modes... As of right now, I am not sure what is causing the odd bump at about 300 μHz.
* **Er_v.gif** : A slightly jarring but demonstrative gif showing the evolution of the radial displacement for frequencies ranging from 40 μHz to 4000 μHz.

## Authors

* **Bridget Andersen**
