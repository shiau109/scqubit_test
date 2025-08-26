import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from qutip import Qobj

omega_onres = 4.5
g = 0.080  # coupling strength

delta_factor = np.linspace(-10, 10, 5)

# Pauli matrices and identity
sx = sigmax()
sy = sigmay()
sz = sigmaz()
I = qeye(2)

V = tensor(sx, sx)


# Arrays to store eigenvalues for plotting
all_evals0 = []
all_evals = []

for df in delta_factor:
    delta = g * df
    omega_0 = omega_onres
    omega_1 = omega_onres + delta
    H0_q0 = Qobj(np.diag([0,omega_0]))
    H0_q1 = Qobj(np.diag([0,omega_1]))
    # Uncoupled Hamiltonian
    H0 = tensor(H0_q0, I) + tensor(I, H0_q1)
    # Full Hamiltonian
    H = H0 + g * V

    # Diagonalize both
    evals0, evecs0 = H0.eigenstates()
    # print(evals0,evecs0)
    print(evals0)

    evals, evecs = H.eigenstates()

    # Store eigenvalues for plotting
    all_evals0.append(evals0)
    all_evals.append(evals)

    # Compute projection (overlap) matrix
    overlap_matrix = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            overlap_matrix[i, j] = abs((evecs0[i].dag() * evecs[j]))**2

    # print(f"\nUncoupled Hamiltonian for delta_factor={df}:")
    # print(H0)
    # print("Full Hamiltonian:")
    # print(H)
    print("Overlap (projection) matrix:")
    print(overlap_matrix)

# Convert lists to numpy arrays for plotting
all_evals0 = np.array(all_evals0)
all_evals = np.array(all_evals)

# Plotting
plt.figure(figsize=(8,6))
for i in range(4):
    plt.plot(delta_factor, all_evals0[:,i], 'o--', label=f'evals0[{i}]')
    plt.plot(delta_factor, all_evals[:,i], 's-', label=f'evals[{i}]')
plt.xlabel('delta_factor')
plt.ylabel('Eigenvalues')
plt.title('Eigenvalues vs delta_factor')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


