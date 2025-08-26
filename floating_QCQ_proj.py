
import scqubits as scq
import numpy as np
import matplotlib.pyplot as plt


test_circuit = scq.Circuit("2FQ1FC.yaml", from_file=True)
# test_circuit.sym_lagrangian()

# test_circuit.sym_hamiltonian()
print(test_circuit.cutoff_names)

print(test_circuit.hamiltonian().shape)
test_circuit.cutoff_n_1 = 10
test_circuit.cutoff_n_2 = 10
test_circuit.cutoff_n_3 = 10
print(test_circuit.hamiltonian().shape)
# print(test_circuit.hamiltonian())
eigenvals = test_circuit.eigenvals()
print(eigenvals)
energy_diff = np.diff(eigenvals)
print(energy_diff)
print(f"E01 {energy_diff[0]}")
print(f"E12 {energy_diff[1]}")

print(f"anharmonicity {(energy_diff[1] -energy_diff[0])}")

Ej = np.linspace(13, 18, 31)

specdata = test_circuit.get_spectrum_vs_paramvals(param_name='qc_EJ1', param_vals=Ej, evals_count=11,subtract_ground=True ) 
specdata.plot_evals_vs_paramvals(subtract_ground=True)
print(specdata.energy_table.shape)

print(specdata.energy_table[0])

plt.show()