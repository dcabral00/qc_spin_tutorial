conda env create -f qc_environment.yml
# Ensure the correct conda environment is active, before using pip
conda activate qc_spin_dynamics_tutorial_env 

# Ensure that the conda version of pip is the one being used
# to enforce that qiskit is contained within that environment
# Then run the following installation commands
# pip install 'qiskit==0.45.0'
# pip install 'qiskit-aer==0.13.0'
# pip install 'qiskit-ibm-provider==0.7.2'

# If ttpy is needed for the classical simulations;
# this fork fixes some issues with the original repository maintained
# by Ivan Oseledets at git+https://github.com/oseledets/ttpy.git 
# pip install git+https://github.com/bcallen95/ttpy.git
