from qiskit.quantum_info import SparsePauliOp


# -----------------------------------------------------------------
# Hamiltonian Functions 
# -----------------------------------------------------------------
def get_hamiltonian_n_site_terms(n, coeff, n_qubits):
    XX_coeff = coeff[0]
    YY_coeff = coeff[1]
    ZZ_coeff = coeff[2]
    Z_coeff = coeff[3]

    XX_term = SparsePauliOp(("I" * n + "XX" + "I" * (n_qubits - 2 - n)))
    XX_term *= XX_coeff
    YY_term = SparsePauliOp(("I" * n + "YY" + "I" * (n_qubits - 2 - n)))
    YY_term *= YY_coeff
    ZZ_term = SparsePauliOp(("I" * n + "ZZ" + "I" * (n_qubits - 2 - n)))
    ZZ_term *= ZZ_coeff
    Z_term = SparsePauliOp(("I" * n + "Z" + "I" * (n_qubits - 1 - n)))
    Z_term *= Z_coeff

    return (XX_term + YY_term + ZZ_term + Z_term)


def get_heisenberg_hamiltonian(n_qubits, coeff=None):

    # Three qubits because for 2 we get H_O = 0
    assert n_qubits >= 3

    if coeff == None:
        coeff = [[1.0, 1.0, 1.0, 1.0] for i in range(n_qubits)]

    # Even terms of the Hamiltonian
    # (summing over individual pair-wise elements)
    H_E = sum((get_hamiltonian_n_site_terms(i, coeff[i], n_qubits)
                for i in range(0, n_qubits-1, 2)))

    # Odd terms of the Hamiltonian
    # (summing over individual pair-wise elements)
    H_O = sum((get_hamiltonian_n_site_terms(i, coeff[i], n_qubits)
                for i in range(1, n_qubits-1, 2)))

    # adding final Z term at the Nth site
    if (n_qubits % 2) == 0:
        H_E += SparsePauliOp(("I" * (n_qubits - 1)  + "Z")) * coeff[n_qubits-1][3]
    else: # SparsePauliOp(("I" * i + "Z" + "I" * (n_qubits - 1 - i))) * coeff[3]
        H_O += SparsePauliOp(("I" * (n_qubits - 1)  + "Z")) * coeff[n_qubits-1][3]

    # Returns the list of the two sets of terms
    return [H_E, H_O]


if __name__ == '__main__':
    num_q = 4
    # XX YY ZZ, Z
    ham_coeffs = ([[0.75/2, 0.75/2, 0.0, 0.65]]
                    + [[0.5, 0.5, 0.0, 1.0]
                    for i in range(num_q-1)])

    spin_chain_hamiltonian = get_heisenberg_hamiltonian(num_q,
                                                        ham_coeffs)
    print('Hamiltonian separated into even and odd components:')
    print(spin_chain_hamiltonian)
    print('Hamiltonian combining even and odd components:')
    print(sum(spin_chain_hamiltonian))
