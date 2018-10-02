# PHYS 643 - Computational Exercise 1

The goal of this computational exercise was to investigate the properties of T = 0 white dwarf stars by integrating the equations of hydrostatic balance.

## Directory contents
* **hw1.pdf** : contains the write-up for this assignment
* **hw1.cc** : the C++ file where the brunt of the RK4 simulation is located
* **pc_mass_radius.txt** : a file containing data mass, radius, density, Fermi energy, and polytropic index data from the simulation
* **mass_radius_relation_plotter.py** : a python script used to plot the mass-radius relation
  * *mass_radius_relation.pdf* : a plot of the mass-radius relation
* **density_profile_plotter.py** : a python script used to plot the density profiles
  * *density_profile_0019.pdf, density_profile_0199.pdf, density_profile_0099.pdf, density_profile_0499.pdf* : a few select plots demonstrating the density profile variation for different masses
  * *density_profile.gif* : a gif demonstrating the density profile evolution for a wide range of masses
* **fermi_energy_plotter.py** : a python script used to plot the Fermi energy vs mass and density
  * *fermi_energy_density.pdf, fermi_energy_mass.pdf* : plots of the Fermi energy vs mass and density
* **polytop_index_plotter.py** : a python script used to plot the polytropic index vs mass and density
  * *gamma_density.pdf, gamma_mass.pdf* : plots of the polytropic index vs mass and density

## Authors

* **Bridget Andersen**

## Acknowledgments

* Much of the derivations were expanded upon the class notes.
* I used [Wikipedia](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge%E2%80%93Kutta_method) to remember how to do Runge Kutta.
