# PHYS 643 - Computational Exercise 2

The goal of this computational exercise was to create a simple 1D hydrodynamic code and use it to explore the properties of a steepening a sound wave.

## Directory contents
* **hw2.pdf** : contains the write-up for this assignment
* **hw2_general.py** : the python file where the donor cell 1D hydrodynamics simulation code is kept
  * Note that the code file [hw2_general.py](https://github.com/bandersen441/PHYS-643/blob/master/Computational_Exercise_2/hw2_general.py) can take command line arguments that specify the particular simulation to run. By default, the code runs a simulation of a low-amplitude gaussian density profile and a zero initial velocity profile. To see options from the command line, type `python hw2_general.py --help`. Feel free to play with the parameters, but I can't promise that the plots will scale nicely by default (I have hardcoded nice scalings for the example simulations). Set `--example_scaling` equal to 0 if you want automatic scaling.

## Authors

* **Bridget Andersen**

## Acknowledgments

* Much of the derivations were expanded upon the class notes.
