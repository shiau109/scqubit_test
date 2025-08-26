import scqubits as scq
 
import numpy as np
import matplotlib.pyplot as plt



circuit_yaml = """# zero-pi
branches:
- [JJ, 1,2, EJ = 20, 0.1fF]
- [C, 1,2, 0.185]
"""


test_circuit = scq.Circuit(circuit_yaml, from_file=False)
test_circuit.sym_lagrangian()

test_circuit.sym_hamiltonian()
# print(test_circuit.hamiltonian())
test_circuit.cutoff_n_1 = 10
eigenvals = test_circuit.eigenvals()
energy_diff = np.diff(eigenvals)
print(energy_diff)
print(f"E01 {energy_diff[0]}")
print(f"E12 {energy_diff[1]}")

print(f"anharmonicity {(energy_diff[1] -energy_diff[0])}")
Ej = np.linspace(10, 30, 20)
# specdata = test_circuit.get_spectrum_vs_paramvals(param_name='EJ', param_vals=Ej ) 
# specdata.plot_evals_vs_paramvals(subtract_ground=True)


transmon = scq.Transmon(EJ=20,
                              EC=0.185,
                              ng=0.0,
                              ncut=10)

# print( transmon.hamiltonian().shape )
print( transmon.hamiltonian() )
print( transmon.eigenvals() )
# print( transmon.eigensys() )
transmon.plot_matrixelements('n_operator', evals_count=20, mode='real')
transmon.plot_matrixelements('n_operator', evals_count=20, mode='imag')

# specdata = transmon.get_spectrum_vs_paramvals(param_name='EJ', param_vals=EJ ) 
# specdata.plot_evals_vs_paramvals()
plt.show()