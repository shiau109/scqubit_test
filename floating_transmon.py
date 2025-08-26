
import scqubits as scq
import numpy as np
import matplotlib.pyplot as plt

circuit_yaml = """
branches:
- [JJ, 1,2, EJ1 = 20, 0.1fF]
- [C, 0,1, Ec_01=97fF]
- [C, 0,2, Ec_02=97fF]
- [C, 1,2, Ec_12=56fF]
"""
file_path = "float_transmon.yaml"
# test_circuit = scq.Circuit(circuit_yaml, from_file=False)
test_circuit = scq.Circuit(file_path, from_file=True)

test_circuit.sym_lagrangian()

test_circuit.sym_hamiltonian()
# print(test_circuit.hamiltonian())
test_circuit.cutoff_n_1 = 20
eigenvals = test_circuit.eigenvals()
print(eigenvals)
energy_diff = np.diff(eigenvals)
print(energy_diff)
print(f"E01 {energy_diff[0]}")
print(f"E12 {energy_diff[1]}")

print(f"anharmonicity {(energy_diff[1] -energy_diff[0])}")

Ej = np.linspace(10, 30, 20)

specdata = test_circuit.get_spectrum_vs_paramvals(param_name='EJ1', param_vals=Ej ) 
specdata.plot_evals_vs_paramvals(subtract_ground=True)


plt.show()