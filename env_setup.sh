conda env create -f qc_environment.yml
# ensure the correct conda environment is active, before using pip
conda activate qc_spin_dynamics_tutorial_env 
# ensure that the conda version of pip is the one being used
# to ensure that qiskit is contained within that environment
# Then run the following installation commands
# pip install 'qiskit==0.45.0'
# pip install 'qiskit-aer==0.13.0'
# pip install 'qiskit-ibm-provider==0.7.2'
