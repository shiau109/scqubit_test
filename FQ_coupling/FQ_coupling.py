import numpy as np
import pandas as pd
import qutip as qt
import matplotlib.pyplot as plt
from scipy.linalg import block_diag
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams['font.family'] = "STIXGeneral"
plt.rcParams['font.size'] = 12


# Calculate coupling strength g_ij between floating transmons qubits from capacitance matrix
# Reference: http://dx.doi.org/10.1103/PhysRevApplied.15.064063 (APPENDIX B)

# Input: Capacitance matrix C (in fF) from Ansys Q3D
# SignalNet: GND, Q0_L (pad1), Q0_R (pad2), Q0_xy (xy line), Q0_read (readout line), C0_L, ......
# where Q* represents qubit, and C* represents coupler


def load_maxwell_C_matrix(filename):
    C = pd.read_csv(filename, skiprows=8, index_col=0)
    N = C.shape[1] - 1
    C = pd.read_csv(filename, skiprows=8, index_col=0, nrows=N)

    # Pre-process to get rid of the stray capacitance to infinity
    for i in range(N):
        C.iloc[i, i] -= np.sum(C.iloc[i, :])
    return C


def calculate_freq_and_gij(C, EJ, print_result=False): # EJ in GHz
    e = 1.60217657e-19  # electron charge
    h = 6.62606957e-34  # Plank's
    hbar = 1.0545718E-34  # Plank's reduced
    # Only keep columns and rows with names ending in '_L' or '_R'
    selected = [name for name in C.columns if name.endswith('_L') or name.endswith('_R')]
    C_matrix = C.loc[selected, selected].to_numpy()

    # Transfrom capacitance matrix to remove the redundant DOF
    num_blocks = C_matrix.shape[0] // 2
    U = block_diag(*([np.array([[1, -1], [1, 1]])]*num_blocks))

    C_reduced = np.linalg.inv(U.T) @ C_matrix @ np.linalg.inv(U)
    qubits = [name.split('_')[0] for name in selected[::2]]

    # Inverse of reduced capacitance matrix
    C_inv = np.linalg.inv(C_reduced)[::2, ::2]
    Ec_matrix = e**2 / (2 * h) * C_inv * 1e6 # Ec in GHz, C in fF
    # print(pd.DataFrame(Ec_matrix * 1e3, index=qubits, columns=qubits))
    
    # Calculate qubit frequencies
    Ec = Ec_matrix.diagonal()
    zeta = (2*Ec / EJ)**0.5
    freq = np.sqrt(8 * EJ * Ec) - Ec * (1 + zeta / 4)


    # Calculate coupling strengths
    g_ij = Ec_matrix / 2**0.5 * ((EJ / Ec)**0.25)[:, None] * ((EJ / Ec)**0.25)[None, :] * (1 - zeta[:, None] / 8 - zeta[None, :] / 8)
    for i in range(len(qubits)):
        g_ij[i, i] = 0
    g_ij *= 2 # (Q1 Q2) * g_ij * (Q1 Q2)^T = g_11 * Q1^2 + g_22 * Q2^2 + 2*g_12*Q1*Q2

    if print_result:
        # print(pd.DataFrame(Ec * 1e3, index=qubits, columns=['Ec (MHz)']))
        print(pd.DataFrame(np.transpose([Ec * 1e3, freq]), index=qubits, columns=['Ec (MHz)', 'Frequency (GHz)']))
        print("\nCoupling strengths g_ij (MHz):")
        print(pd.DataFrame(g_ij * 1e3, index=qubits, columns=qubits))

    return freq, g_ij, Ec, zeta

def get_eig(C, EJ):
    freq, g_ij, Ec, zeta = calculate_freq_and_gij(C, EJ)

    dim = 10  # number of levels per qubit
    N = np.size(freq)
    H = 0

    # \sum_i \omega_i (a^dagger a) + Ec_i/2 (1 + zeta_i/4 - (1 + 9*zeta_i/16) a^dagger a) a^dagger a
    for i in range(N):
        H_sub = (freq[i] + Ec[i] / 2 * ((1 + zeta[i] / 4) - (1 + 9 * zeta[i] / 16) * qt.num(dim))) * qt.num(dim)
        H += qt.tensor([H_sub if j == i else qt.qeye(dim) for j in range(N)])
    
    # \sum_{i<j} -g_ij (a^dagger - a)(b^dagger - b)
    for i in range(N): 
        for j in range(i+1, N):
            H -= g_ij[i, j] * qt.tensor([ (qt.create(dim) - qt.destroy(dim)) if k == i else (qt.create(dim) - qt.destroy(dim)) if k == j else qt.qeye(dim) for k in range(N)])
    
    eigenvalues = H.eigenenergies()
    return eigenvalues

def calculate_readout(C, EJ, fr):
    e = 1.60217657e-19  # electron charge
    h = 6.62606957e-34  # Planck's
    hbar = 1.0545718E-34  # Planck's reduced

    names = [name for name in C.columns if name.endswith('_L')]
    qubits = [name.split('_')[0] for name in names]
    for q in qubits:
        if f"{q}_read" not in C.columns:
            continue
        selected = [f"{q}_read", f"{q}_L", f"{q}_R"]
        Cr = 0.25 * np.pi / (2*np.pi * fr) / 50 * 1e6 #For lambda/4 resonator (fF)
        # Cr = 0.5 * np.pi / (2*np.pi * fr) / 50 * 1e6 #For lambda/2 resonator (fF)

        C_matrix = C.loc[selected, selected].to_numpy()
        C_matrix[0, 0] += C.loc[f"{q}_read", "GND"] + Cr # Add resonator capacitance to ground

        U = block_diag(*([1] + [np.array([[1, -1], [1, 1]])]))
        C_transform = np.linalg.inv(U.T) @ C_matrix @ np.linalg.inv(U)
        C_inv = np.linalg.inv(C_transform)[:2, :2]
        Ec_matrix = e**2 / (2 * h) * C_inv * 1e6 # Ec in GHz, C in fF

        # Calculate qubit frequencies
        Ec = Ec_matrix[1, 1]
        zeta = (2*Ec / EJ)**0.5
        freq = np.sqrt(8 * EJ * Ec) - Ec * (1 + zeta / 4)

        # Calculate coupling strengths
        g = Ec_matrix[0, 1] / 2**0.5 * (EJ / Ec)**0.25 * (fr**2 / 2 / (e**2 / Cr / h * 1e6)**2)**0.25 * (1 - zeta / 8)
        g *= 2 # (Q1 Q2) * g_ij * (Q1 Q2)^T = g_11 * Q1^2 + g_22 * Q2^2 + 2*g_12*Q1*Q2
        print(f"Qubit {q}: freq = {freq:.3f} GHz, Ec = {Ec*1e3:.1f} MHz, g = {g*1e3:.1f} MHz")
        

if __name__ == '__main__':
    filename = "Capacitance_2x2FQ.csv"
    #filename = r"C:\Users\linyc\Documents\Q3D\2FQ1FC-H.csv"

    C = load_maxwell_C_matrix(filename)
    freq, g_ij, Ec, zeta = calculate_freq_and_gij(C, EJ=16, print_result=True)
    calculate_readout(C, 14, 6)

    EJ_c0 = np.linspace(8, 16, 81)  # in GHz
    Eig = []
    for i in range(len(EJ_c0)):
        EJ = np.array([EJ_c0[i], 10, 14])
        Eig.append(get_eig(C, EJ))
    Eig = np.array(Eig)
    plt.figure()
    for i in range(10):
        plt.plot(EJ_c0, Eig[:, i], label=f'Level {i}')
    plt.xlabel('EJ_c0 (GHz)')
    plt.ylabel('Eigenenergy (GHz)')
    plt.title('Eigenenergies vs EJ')
    plt.tight_layout()
    print(np.min(Eig[:, 2] - Eig[:, 1]))
    plt.show()
