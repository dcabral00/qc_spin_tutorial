The purpose of this repository is to give a barebones introduction
to performing quantum dynamics simulations using a qubit-based
platform. 

Most of the functions are documented and grouped based on their
usage. 

circuit_preliminaries.py - walks through some of the initial steps 
    on how to set up quantum circuits, printing them, and perform
    measurements. It also demonstrates the circuit used for time
    evolution
hamiltonians.py - contains hamiltonian used in the simulations
propagators.py - contains the dynamic propagtors uses to evolve
    the initial state

The main simulation modules are listed below; the output of executing
them is both plotted as well as saved for user convenience.

statevector_simulation.py - contains the function responsible for
    performing statevector simulation (qsolve) by reinitialization
    and when executed generates the plot for a given set of
    parameters set in the figure.
hadamard_test.py - contains the functions responsible for executing
    the hadamard test and post process the results. (takes about 1
    hour for the parameters included by default using Intel i7 12k
    processor when run locally)

For the sake of observing the results of execution without running
the simulation code the plot_results module can be used. By default
pre-generated data is provided within the data/ folder and loaded
as needed.

plot_results.py - plots the pre-generated files

Extra files include plot_style.txt to make pretty python figures.
The default color pallete is colorblind friendly.

Disclaimer: While it is possible to run these on a jupyter notebook,
the author recommends using the included python modules, in 
particular if the user wants to take the next step in extending them
to make production ready code for research (as the author has done). 
The author hope that this serves as a template in doing so, which was
not easily available when starting research in the field.

Furthermore, the author welcomes feedback and enhancement of this
repository to make it more readable and accessible for pedagogical
purposes.

