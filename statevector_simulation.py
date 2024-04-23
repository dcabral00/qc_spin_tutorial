import numpy as np
import matplotlib.pyplot as plt
from qiskit import execute
from qiskit import QuantumCircuit, QuantumRegister
from qiskit import BasicAer
from propagators import get_time_evolution_operator


dpi=300
plt.style.use('plot_style.txt')
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['lines.markersize'] = 11 
plt.rcParams["figure.figsize"] = (6.4, 3.6)


# Quantum circuit for propagation
def qsolve_statevector(psin, qc):
    '''
        Performs dynamical propagation for a specific quantum circuit, using the
        an initial state equal to the result of the prior iteration. In effect,
        this propagates:

        | \psi _t \rangle  = e^{i*\tau*H/hbar} e^{i*\tau*H/hbar} ... | \psi _0 \rangle
        -> | \psi _t \rangle  = e^{i*\tau*H/hbar} | \psi _{t-\tau} \rangle

        Use of reinitialization in this way is possible due to the usage of the
        statevector_simulator backend; this may not be practical for an
        implementation on actual hardware.
        Returns:
            psin (statevector): final statevector after execution
    '''
    # Determining number of qubits from the length of the state vector
    d=int(np.log2(np.size(psin)))
    # Circuit preparation
    qre = QuantumRegister(d)
    circ = QuantumCircuit(qre)
    circ.initialize(psin,qre)
    circ.barrier()
    circ.append(qc, qre)
    circ.barrier()
    # Circuit execution
    device = BasicAer.get_backend('statevector_simulator')
    psin = execute(circ, backend=device).result()
    return psin.get_statevector()


if __name__ == '__main__':

    num_q = 3
    evolution_timestep = 0.1
    n_trotter_steps = 1
    # XX YY ZZ, Z
    hamiltonian_coefficients = ([[0.75/2, 0.75/2, 0.0, 0.65]]
                                + [[0.5, 0.5, 0.0, 1.0]
                                    for i in range(num_q-1)])
    num_shots = 100

    # Qubit basis states
    zero_state = np.array([[1],[0]])
    one_state = np.array([[0],[1]])

    # For a 011 initial state prepare as follows
    psin = zero_state # for the first spin
    # iterates over the remaining spins, by performing
    # Kronecker Product
    for i in range(num_q-1):
        psin = np.kron(psin, one_state)
    psin0 = psin.flatten()
    print(psin0)

    # time evolution operator
    time_evo_op = get_time_evolution_operator(num_qubits=num_q,
            tau=evolution_timestep,
            trotter_steps=n_trotter_steps,
            coeff=hamiltonian_coefficients)
    # number of steps for which to propagate
    # (totaling 25 units of time)
    nsteps = 250
    psin_list = []
    psin_list.append(psin0)
    correlation_list = []
    # performs dynamical propagation by statevector re-initialization
    for k in range(nsteps):
        print(f'Running dynamics step {k}')
        if k > 0:
            psin = qsolve_statevector(psin_list[-1], time_evo_op)
            # removes the last initial state to save memory
            psin_list.pop()
            # stores the new initial state
            psin_list.append(psin)
        correlation_list.append(np.vdot(psin_list[-1],psin0))
    
    time = np.arange(0, evolution_timestep*(nsteps),
                     evolution_timestep)
    np.save(f'{num_q}_spin_chain_time', time)
    sa_observable = np.abs(correlation_list)
    np.save(f'{num_q}_spin_chain_SA_obs', sa_observable)
    # plotting
    plt.plot(time, sa_observable, '-o')
    plt.xlabel('Time')
    plt.ylabel('Absolute Value of Survival Amplitude, '
               r'$\left|\langle \psi | \psi \rangle \right|$')
    plt.xlim((min(time), max(time)))
    plt.yscale('log')
    plt.legend()
    plt.savefig('statevector.pdf', format='pdf', dpi=dpi)
    plt.savefig('statevector.png', format='png', dpi=dpi)
    plt.show()
