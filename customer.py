
import scqubits as scq
import numpy as np
import matplotlib.pyplot as plt

circuit_yaml = """# zero-pi
branches:
- [JJ, 1,2, EJ1 = 20, 0.1fF]
- [C, 0,1, Ec_01=0.199]
- [C, 0,2, Ec_02=0.199]
- [C, 1,2, Ec_12=0.343]
"""
print(1/(1/0.343+1/(0.199*2)))

test_circuit = scq.Circuit(circuit_yaml, from_file=False)
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