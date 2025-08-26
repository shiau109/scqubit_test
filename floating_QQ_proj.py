
import scqubits as scq
import numpy as np
import matplotlib.pyplot as plt

circuit_yaml = """# zero-pi
branches:
- [JJ, 0,1, q0_EJ1 = 20, 0.01fF]
- [C, 0,1, q0_Ec_01=100fF]
- [JJ, 0,2, q1_EJ1 = 15, 0.01fF]
- [C, 0,2, q1_Ec_01=100fF]
- [C, 1,2, c01=5fF]

"""

test_circuit = scq.Circuit(circuit_yaml, from_file=False)

test_circuit.sym_lagrangian()
test_circuit.sym_hamiltonian()
test_circuit.sym_lagrangian(vars_type="new")

print(test_circuit.cutoff_names)
print(test_circuit.var_categories)


test_circuit.cutoff_n_1 = 10
test_circuit.cutoff_n_2 = 10
# print(test_circuit.hamiltonian())
import seaborn as sns

# Imaginary part
# plt.figure()
# sns.heatmap(test_circuit.hamiltonian().toarray().imag, annot=False, cmap="coolwarm", cbar=True)
# plt.title("csc_matrix imaginary values")

# Real part
plt.figure()
sns.heatmap(test_circuit.hamiltonian().toarray().real, annot=False, cmap="coolwarm", cbar=True)
plt.title("csc_matrix real values")

# test_circuit.plot_matrixelements('hamiltonian', evals_count=49)

# print(test_circuit.branches)
# print(test_circuit.transformation_matrix)



# system_hierarchy = [[0,1], [2,3]]
# print(test_circuit.hamiltonian().shape)
# print(scq.truncation_template(system_hierarchy))
# test_circuit.configure(system_hierarchy=system_hierarchy)
# print(test_circuit.hamiltonian().shape)

# # print(test_circuit.hamiltonian())
# eigenvals = test_circuit.eigenvals()
# print(eigenvals)
# energy_diff = np.diff(eigenvals)
# print(energy_diff)
# print(f"E01 {energy_diff[0]}")
# print(f"E12 {energy_diff[1]}")

# print(f"anharmonicity {(energy_diff[1] -energy_diff[0])}")

Ej = np.linspace(15, 25, 21)
g01 = np.linspace(1, 3, 11)
# specdata = test_circuit.get_spectrum_vs_paramvals(param_name='g01', param_vals=g01, evals_count=6 ) 
specdata = test_circuit.get_spectrum_vs_paramvals(param_name='q1_EJ1', param_vals=Ej, evals_count=9 ) 

specdata.plot_evals_vs_paramvals(subtract_ground=True)


plt.show()