
import scqubits as scq
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr




# --- Function to save energy table using xarray ---
def save_energy_table_xr(energy_table, Ej, filename):
	ds = xr.Dataset({
		'energy': (['Ej', 'level'], energy_table)
	}, coords={'Ej': Ej, 'level': np.arange(energy_table.shape[1])})
	ds.to_netcdf(filename)
	print(f"Energy table saved to {filename} (NetCDF/xarray)")

# --- Function to load and plot energy table from xarray ---
def plot_energy_table_from_xr(filename):
    ds = xr.load_dataset(filename)
    # Calculate chi
    spectra_T = ds.energy.transpose()
    chi = spectra_T.sel(level=8) - (spectra_T.sel(level=3) + spectra_T.sel(level=4))
    print(chi)

    # Plot chi vs Ej
    plt.figure(figsize=(8,6))
    plt.plot(ds.Ej.values, chi, 'o-', label='chi')
    plt.xlabel('qc_EJ1')
    plt.ylabel('Chi')
    plt.title('Chi vs qc_EJ1')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Plot all energy levels
    plt.figure(figsize=(8,6))
    for i in ds.level.values:
        plt.plot(ds.Ej.values, ds.energy.sel(level=i), label=f'Level {i}')
    plt.xlabel('qc_EJ1')
    plt.ylabel('Energy')
    plt.title('Loaded Energy Table Spectrum (xarray)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


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

Ej = np.linspace(5, 15, 101)

specdata = test_circuit.get_spectrum_vs_paramvals(param_name='qc_EJ1', param_vals=Ej, evals_count=11,subtract_ground=True ) 
print(specdata.energy_table.shape)


save_energy_table_xr(specdata.energy_table, Ej, "energy_table.nc")


# Example usage: plot from saved file
plot_energy_table_from_xr("energy_table.nc")

plt.show()