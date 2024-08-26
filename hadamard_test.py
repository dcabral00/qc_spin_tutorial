import numpy as np
import matplotlib.pyplot as plt
from custom_execute import execute
from qiskit_aer import Aer
from qiskit import QuantumCircuit
from qiskit import QuantumRegister, ClassicalRegister
from propagators import get_time_evolution_operator


def get_hadamard_test(num_q, initial_state, control_operation,
                      control_repeats=0, imag_expectation=False):
    '''
        Writes the Hadamard test circuit to evaluate the real (or imaginary)
        components of a given operator (or set of operators) by measuring the
        ancilla.

        Inputs:
            num_qubits (int): number of qubits needed for the calculation
            initial_state: quantum circuit object containing the initialization
                for the system
            control_operation: quantum circuit object containing the controlled
                version of the unitary operator for which to evaluate the
                expectation value
            control_repeats (int): number of repetitions of the controlled
                operation in the Hadamard test circuit (such as time-propagation)
            imag_expectation (boolean): whether to assemble the circuit for
                evaluation of the real or imaginary components of the operator.
        Returns:
            qc_hadamard: quantum circuit object containing the hadamard test
    '''
    # Circuit object framework
    qr_hadamard = QuantumRegister(num_q+1)
    cr_hadamard = ClassicalRegister(1)
    qc_hadamard = QuantumCircuit(qr_hadamard, cr_hadamard) # instantiated here
    # Initialization of calculation qubits
    qc_hadamard.append(initial_state, qr_hadamard[1:]) # initial psi
    qc_hadamard.barrier()
    # Hadamard test structure
    qc_hadamard.h(0)
    if imag_expectation:
        qc_hadamard.p(-np.pi/2, 0) # qc_hadamard.s(0).inverse() may be equivalent
    # iterates over the number of times the control operation should be added
    for i in range(control_repeats):
        qc_hadamard.append(control_operation, qr_hadamard[:])
    qc_hadamard.h(0)
    qc_hadamard.barrier()
    # Measuring the ancilla
    qc_hadamard.measure(0,0)
    return qc_hadamard


def get_circuit_execution_counts(qc, backend, n_shots=100):
    '''
        Takes a quantum circuit, Qiskit supported backend and a specified number
        of execution times to calculate the number of times each state in the
        circuit was measured.

        Inputs:
            qc: Qiskit quantum circuit object
            backend: qiskit supported quantum computer framework (either simulator
                or actual device)
                n_shots (int): default=100; number of times for which to execute
                    the circuit using the specified backend.
    '''
    qc_execution = execute(qc, backend, shots=n_shots)
    counts = qc_execution.result().get_counts()
    return counts # number of times measured 0 and 1


def get_spin_correlation(counts):
    '''
        Generates the spin correlation at a particular time by averaging over
        the counts of the ancillary qubit.
    '''

    qubit_to_spin_map = {
        '0': 1,
        '1': -1,
    }
    total_counts = 0
    values_list = []
    for k,v in counts.items():
        values_list.append(qubit_to_spin_map[k] * v)
        total_counts += v
    # print(values_list)
    average_spin = (sum(values_list)) / total_counts
    return average_spin


def get_initialization(num_qubits, initialization_string):
    '''
        Creates a circuit containing the qubit initialization for a specified
        number of qubits and given initialization string. There are alternative
        methods to initialize the circuit besides the one used in this function
        including using a statevector or amplitude embedding of an arbitrary state.
    '''
    qr_init = QuantumRegister(num_qubits)
    qc_init = QuantumCircuit(qr_init)
    qc_init.initialize(initialization_string, qr_init[:]) # 0x1x1
    return qc_init


if __name__ == '__main__':
    # switch this variable to True to execution the 'for' loop.
    # it takes >1hr for 3 spins, with the parameters defined above
    # and 1000 shots for each time step.
    # run_hadamard_test_boolean = False
    run_hadamard_test_boolean = True

    # IMPORTANT: Use qasm_simulator to obtain meaningful statistics
    # statevector is not appropriate for this method
    simulator = Aer.get_backend('qasm_simulator')

    num_q = 3
    n_trotter_steps = 1
    # XX YY ZZ, Z
    hamiltonian_coefficients = ([[0.75/2, 0.75/2, 0.0, 0.65]]
                                + [[0.5, 0.5, 0.0, 1.0]
                                    for i in range(num_q-1)])

    num_shots = 100 # increase to check for convergence

    evolution_timestep=0.1
    total_time = 250
    time_range = np.arange(0, total_time+evolution_timestep,
                           evolution_timestep)


    # time evolution operator
    time_evo_op = get_time_evolution_operator(num_qubits=num_q,
            tau=evolution_timestep,
            trotter_steps=n_trotter_steps,
            coeff=hamiltonian_coefficients)
    # for closed controlled operation (ie if qubit is 1)
    controlled_time_evo_op = time_evo_op.control()
    # for open controlled operation (ie if qubit is 0);
    # not used here for this correlation fxn but needed for others
    # controlled_time_evo_op = time_evo_op.control(ctrl_state=0)
    print(controlled_time_evo_op.decompose())


    init_state_list = '1' + '0' * (num_q-1)
    init_circ = get_initialization(num_q, init_state_list)
    init_circ.draw(style='iqp')
    print(init_circ)

    # lists t store observables
    if run_hadamard_test_boolean:
        real_amp_list = []
        imag_amp_list = []
        for idx,time in enumerate(time_range):
            print(f'Running dynamics step {idx}')
            # Real component ------------------------------  
            qc_had_real = get_hadamard_test(num_q, init_circ,
                                            controlled_time_evo_op,
                                            control_repeats=idx,
                                            imag_expectation=False)
            had_real_counts = get_circuit_execution_counts(
                    qc_had_real, simulator, n_shots=num_shots)
            real_amplitude = get_spin_correlation(had_real_counts)
            real_amp_list.append(real_amplitude)

            # Imag component ------------------------------  
            qc_had_imag = get_hadamard_test(num_q, init_circ,
                                            controlled_time_evo_op,
                                            control_repeats=idx,
                                            imag_expectation=True)
            had_imag_counts = get_circuit_execution_counts(
                    qc_had_imag, simulator, n_shots=num_shots)
            imag_amplitude = get_spin_correlation(had_imag_counts)
            imag_amp_list.append(imag_amplitude)
            print(f'Finished step {idx}, where '
                  f'Re = {real_amplitude:.3f} '
                  f'Im = {imag_amplitude:.3f}')
    
        real_amp_array = np.array(real_amp_list)
        imag_amp_array = np.array(imag_amp_list)

        print('Saving data...')
        # real and imaginary amplitude components
        np.savetxt('real_amp_array.csv', real_amp_array,
                   fmt='%.18e', delimiter=';', newline='\n')
        np.savetxt('imag_amp_array.csv', imag_amp_array,
                   fmt='%.18e', delimiter=';', newline='\n')
        # survival amplitude
        np.savetxt('np_abs_correlation_with_hadamard_test.csv',
                   np.abs(real_amp_array + 1j*imag_amp_array),
                   fmt='%.18e', delimiter=';', newline='\n')
        # survival probability
        np.savetxt('np_sqrt_sum_squares_correlation_with_hadamard'
                   '_test.csv',
                   np.sqrt(real_amp_array**2 + imag_amp_array**2),
                   fmt='%.18e', delimiter=';', newline='\n')
    real_amp_array = np.loadtxt('real_amp_array.csv',
                                delimiter=';')
    imag_amp_array = np.loadtxt('imag_amp_array.csv',
                                delimiter=';')
    np_abs_correlation_with_hadamard_test = np.loadtxt(
            'np_abs_correlation_with_hadamard_test.csv',
            delimiter=';')
    np_sqrt_sum_squares_correlation_with_hadamard_test = np.loadtxt(
            'np_sqrt_sum_squares_correlation_with_hadamard_test.csv',
            delimiter=';')


    # plotting the data
    plt.plot(time_range, np_abs_correlation_with_hadamard_test,
             '.', label='Hadamard Test')

    sa_statevector = np.load(f'data/{num_q}_spin_chain_SA_obs.npy')
    time = np.load(f'data/{num_q}_spin_chain_time.npy')
    plt.plot(time, sa_statevector, '-', label='Statevector')

    plt.xlabel('Time')
    plt.ylabel('Absolute Value of Survival Amplitude')
    plt.legend()
    plt.show()
