# from qiskit.QuantumInfo import SparsePauliOp
import numpy as np
# from qiskit_aer import Aer
 
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.quantum_info.operators import Operator, Pauli
from qiskit.quantum_info import process_fidelity
 
from qiskit.circuit.library import RXGate, XGate, CXGate


# see https://docs.quantum.ibm.com/build/operators-overview
# note that (I x XX x I) = IxXxI * IxXxI

# def dissipation_1(n_qubits):
#     for n in range(nspins):
#         # identities that may precede site n
#         zn = "I"*n + "Z" + "I"*(nspins-n-1)
#         H_diag += SparsePauliOp((zn)) * h_list[n] * pref
  
  # for n in range(nspins):
  #     # identities that may precede site n
  #     zn = "I"*n + "Z" + "I"*(nspins-n-1)
  #     H_diag += SparsePauliOp((zn)) * h_list[n] * pref
  # return H_diag
    


# # ----------------------------------------------------------------
# # Hamiltonian Routines
# 
# def get_XXZ_qiskit(n_qubits, h_list, J_list, prefactor=True):
#     H_diag = get_diag_XXZ_qiskit(n_qubits, h_list, prefactor)
#     H_off_diag = get_off_diag_XXZ_qiskit(n_qubits, J_list, prefactor)
#     return H_diag + H_off_diag
# 
# 
# def get_diag_XXZ_qiskit(nspins, h_list, prefactor=True):
#     if prefactor:
#         pref = 0.5
#     else:
#         pref = 1
#     H_diag = 0
#     for n in range(nspins):
#         # identities that may precede site n
#         zn = "I"*n + "Z" + "I"*(nspins-n-1)
#         H_diag += SparsePauliOp((zn)) * h_list[n] * pref
#     return H_diag
# 
# 
# def get_off_diag_XXZ_qiskit(nspins, J_list, prefactor=True):
#     if prefactor:
#         pref = 0.5
#     else:
#         pref = 1
#     H_off = 0
#     for n in range(1, nspins):
#         # identities that may follow site n
#         xn = "I"*(n-1) + "XX" + "I"*(nspins-n-1)
#         yn = "I"*(n-1) + "YY" + "I"*(nspins-n-1)
#         zn = "I"*(n-1) + "ZZ" + "I"*(nspins-n-1)
#         H_off += (SparsePauliOp((xn)) * J_list[n] * pref**2
#                   + SparsePauliOp((yn)) * J_list[n] * pref**2
#                   + SparsePauliOp((zn)) * J_list[n] * pref**2)
#         # adding to off-diagonal hamiltonian
#         # H_off += hxn + hyn + hzn
#     return H_off


if __name__ == '__main__':
    print()
