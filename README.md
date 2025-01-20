[![License: GNU AGPL v3](https://img.shields.io/badge/License-GNU_AGPL_v3-lightgrey.svg)](LICENCE)
[![Static Badge](https://img.shields.io/badge/CQDMQD-00268d?style=flat&logoColor=00268d&label=NSF&labelColor=00268d&color=00268d&link=https%3A%2F%2Fcqdmqd.yale.edu%2F)](https://cqdmqd.yale.edu/)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/dcabral00/qc_spin_tutorial/blob/main/qc_spin_chain_tutorial.ipynb)


# Quantum Dynamics Simulations with Qubit-Based Platforms

This repository provides a minimal introduction to performing quantum dynamics simulations using qubit-based platforms. The provided modules are organized by functionality and designed for ease of extension.

## Table of Contents
1. [Getting Started](#getting-started)
   - [Installation](#installation)
     - [Using Conda](#using-conda)
     - [Manual Installation](#manual-installation)
   - [Modules Overview](#modules-overview)
3. [Jupyter Notebook Demo](#jupyter)
4. [Results Visualization](#visualization)
5. [Disclaimer](#disclaimer)
6. [Contact](#contact)
7. [License](#license)

## Getting Started <a name="getting-started"></a>


### Installation <a name="installation"></a>

To set up the environment for executing the code in this repository, you have two options: using `conda` or installing packages manually.

#### Using Conda <a name="using-conda"></a>

1. **Clone the Repository**:
   Choose a target folder location and clone the repository:
   ```bash
   git clone git@github.com:dcabral00/qc_spin_tutorial.git
   cd qc_spin_tutorial
   ```

2. **Create and Activate the Conda Environment**:
   Set up the environment using the provided `.yml` file:
   ```bash
   conda env create -f qc_environment.yml
   conda activate qc_spin_dynamics_tutorial_env
   ```

   This will install all necessary dependencies for running the scripts in this repository.


#### Manual Installation <a name="manual-installation"></a>

Alternatively, you can manually install the required packages listed in the `requirements.txt` file:

1. **Clone the Repository**:
   As above, clone the repository and navigate to the project directory:
   ```bash
   git clone git@github.com:dcabral00/qc_spin_tutorial.git
   cd qc_spin_tutorial
   ```

2. **Install Dependencies**:
   Use `pip` to install the packages:
   ```bash
   pip install -r requirements.txt
   ```

   This command will install all the dependencies specified in the `requirements.txt` file.


### Modules Overview <a name="modules-overview"></a>
#### Preliminary Setup <a name="setup"></a>

`circuit_preliminaries.py` - Walks through setting up quantum circuits, printing, and performing measurements. Demonstrates time-evolution circuits.

`hamiltonians.py` - Contains Hamiltonians used in simulations.

`propagators.py` - Implements dynamic propagators for evolving the initial state.
    
#### Core Modules <a name="core"></a>

`statevector_simulation.py` - Performs statevector simulations (qsolve) by reinitializing. Outputs plots for given parameters.

`hadamard_test.py` - Realistic simulation protocol for Hamiltonian evolution and calculation of observables by executing the Hadamard test.
    
#### Results Visualization <a name="visualization"></a>

`plot_results.py` - Plots pre-generated data stored in the `data/` folder. Useful for visualizing results without running simulations.


## Jupyter Notebook Demo <a name="jupyter"></a>

For a quick and interactive way to explore the usage of this code, use
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/dcabral00/qc_spin_tutorial/blob/main/qc_spin_chain_tutorial.ipynb)

```bash
jupyter notebook qc_spin_chain_tutorial.ipynb
```

## Disclaimer <a name="disclaimer"></a>

While the code can be executed within a Jupyter notebook, we recommend using the provided Python modules for research applications and production-level development. The Jupyter notebook version, while convenient for initial exploration and experimentation, may lack the robustness and scalability of the provided Python modules. A version of this repository has been integrated on a forthcoming publication, whose citation will be included in this repository. This codebase has been used as the basis for research work in preparation and submitted for publication.


## Contact <a name="contact"></a>

Feedback and enhancements to make the repository more accessible for pedagogical purposes are welcome. For questions, comments, or support, please contact:

Delmar G. A. Cabral (delmar.azevedocabral@yale.edu)

## License <a name="license"></a>

This repository is licensed under the GNU AGPL v3. See the `LICENSE` file for details.
