# from qiskit.QuantumInfo import SparsePauliOp
import numpy as np
# from qiskit_aer import Aer
 
from qiskit.quantum_info.operators import Operator, Pauli

# see https://docs.quantum.ibm.com/build/operators-overview
# ----------------------------------------------------------------
# Operators

# Pauli SU(2) matrices
sigmaz = np.array([[1.0, 0], [0, -1]])
sigmay = np.array([[0, -1.0j], [1j, 0]])
sigmax = np.array([[0.0, 1], [1, 0]])
identm = np.array([[1.0, 0], [0, 1.0]])

# Ladder operators ()
sigmam = (sigmax+1j*sigmay)*0.5
sigmap = (sigmax-1j*sigmay)*0.5
sigman = sigmap @ sigmam

def sigmaX():
    return Operator(Pauli('X'))

def sigmaY():
    return Operator(Pauli('Y'))

def sigmaZ():
    return Operator(Pauli('Z'))

def sigmaI():
    return Operator(Pauli('I'))

def sigmaM():
    return 0.5 * (sigmaX() + 1j * sigmaY())

def sigmaP():
    return 0.5 * (sigmaX() - 1j * sigmaY())

def sigmaN():
    return sigmaM().compose(sigmaP())


if __name__ == '__main__':
    print(sigmaX())
    print(sigmaY())
    print(sigmaZ())
    print(sigmaI())
    print(sigmaX())
    print('------------------------------------------------------')
    print(sigmaM())
    print(sigmam)
    print('------------------------------------------------------')
    print(sigmaP())
    print(sigmap)
    print('------------------------------------------------------')
    print(sigmaN())
    print(sigman)
    print('------------------------------------------------------')
