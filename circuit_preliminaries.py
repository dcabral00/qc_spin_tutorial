# used to compile the circuit, tranpiling to the set of
# operations supported by the device backend and execution
from qiskit import QuantumCircuit
from qiskit import QuantumRegister, ClassicalRegister
from qiskit import transpile
from custom_execute import execute 
from qiskit_aer import Aer

from propagators import get_time_evolution_operator


num_q = 3
evolution_timestep = 0.1
n_trotter_steps = 1
# XX YY ZZ, Z
hamiltonian_coefficients = ([[0.75/2, 0.75/2, 0.0, 0.65]]
                            + [[0.5, 0.5, 0.0, 1.0]
                                for i in range(num_q-1)])
num_shots = 100


if __name__ == '__main__':
    # specifying a quantum register with specific number of qubits
    qr = QuantumRegister(num_q)
    # classical register used for measurement of qubits
    cr = ClassicalRegister(num_q)
    # quantum circuit combining quantum and classical registers
    qc = QuantumCircuit(qr, cr) # instantiated here
    qc.draw(style='iqp')
    print(qc)
    
    # specifying initial state by flipping qubit states
    for qubit_idx in range(num_q):
        if qubit_idx == 0:
            # generate only one spin-up at first qubit
            qc.id(qubit_idx)
        else:
            # make all other spins have the spin-down state
            qc.x(qubit_idx)
    qc.barrier()
    qc.draw(style='iqp')
    print(qc)
    
    # checking the initial state
    device = Aer.get_backend('statevector_simulator')
    qc_init_state = execute(qc, backend=device).result()
    qc_init_state = qc_init_state.get_statevector()
    print(qc_init_state)
    
    # generating the time evolution operator for a specific set of
    # hamiltonian parameters and timestep
    time_evo_op = get_time_evolution_operator(num_qubits=num_q,
            tau=evolution_timestep,
            trotter_steps=n_trotter_steps,
            coeff=hamiltonian_coefficients)
    
    # appending the Hamiltonian evolution to the circuit
    qc.append(time_evo_op, list(range(num_q)))
    qc.barrier()
    qc.draw(style='iqp')
    print(qc)

    # Depth check
    print('Depth of the circuit is', qc.depth())
    # transpiled circuit to statevector simulator
    qct = transpile(qc, device, optimization_level=2)
    qct.decompose().decompose()
    qct.draw(style='iqp')
    print(qct)

    print('Depth of the circuit after transpilation is '
          f'{qct.depth()}')
    
    # showcasing simple execution
    qct_run = execute(qct, device, shots=1270).result()
    qct_run_counts = qct_run.get_counts()
    print(f'Circuit measurement counts {qct_run_counts}')
    
