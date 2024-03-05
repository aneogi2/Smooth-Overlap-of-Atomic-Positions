#Import Modules
import os
from glob import glob
from natsort import natsorted
import matplotlib.pyplot as plt
from pymatgen.core.periodic_table import Element, DummySpecie
import numpy as np
import pandas as pd
from itertools import combinations
from dscribe.descriptors import SOAP
from ase.io import read, write
from ase.atoms import Atoms
#from ase.io.lammps import read_lammps_data
from dscribe.core.system import System
from mendeleev import element
from scipy.spatial.distance import pdist, squareform
from scipy.spatial.distance import euclidean

#Read periodic boundary condition lammps files
def read_lammps_data(lammps_data_file):
    with open(lammps_data_file, 'r') as f:
        lines = f.readlines()

    num_atoms = int(lines[1].split()[0])  # Extract the number of atoms directly from the second line

    # Find the line with 'Atoms' keyword to start reading atom data
    atom_data_start = None
    for idx, line in enumerate(lines):
        if line.strip().startswith('Atoms'):  # Adjusted to find 'Atoms' section
            atom_data_start = idx + 2  # Skip the 'Atoms' line and the blank line after
            break

    # Read atom data and store in the 'data' list, ignoring any additional columns beyond the first 5
    data = []
    for line in lines[atom_data_start: atom_data_start + num_atoms]:
        atom_col = line.split()
        atom_type = int(atom_col[1])  # Atom type as integer
        position = [float(val) for val in atom_col[2:5]]  # x, y, z coordinates
        data.append((atom_type, position))

    # Extract PBC from the LAMMPS data file
    xlo, xhi = [float(val) for val in lines[4].split()[0:2]]
    ylo, yhi = [float(val) for val in lines[5].split()[0:2]]
    zlo, zhi = [float(val) for val in lines[6].split()[0:2]]

    # Assuming all atoms are of the same type for simplicity; adjust as necessary
    symbols = ['Si'] * num_atoms  # 'Si' for Slicon; adjust based on actual atom types in your file
    positions = [position for _, position in data]

    atoms_bulk = Atoms(symbols=symbols, positions=positions, pbc=[1, 1, 1])
    atoms_bulk.cell = [[xhi - xlo, 0, 0], [0, yhi - ylo, 0], [0, 0, zhi - zlo]]

    return atoms_bulk

"""Map chemical symbols to a specified atomic number for an Atoms object."""
def map_atomic_numbers(atoms, atomic_number):
    atomic_numbers = [atomic_number] * len(atoms.get_chemical_symbols())
    atoms_with_atomic_numbers = Atoms(symbols=atoms.get_chemical_symbols(),
                                      positions=atoms.get_positions(),
                                      pbc=atoms.get_pbc(),
                                      cell=atoms.get_cell())
    atoms_with_atomic_numbers.set_atomic_numbers(atomic_numbers)
    return atoms_with_atomic_numbers

# Specify the SOAP descriptor parameters
element_numbers = [8]
rcut = 8.0
nmax = 8
lmax = 6

# Create the average SOAP descriptor with "inner" averaging
average_soap = SOAP(
    species=element_numbers,
    r_cut=rcut,
    n_max=nmax,
    l_max=lmax,
    average="inner",
    sparse=False
)

#Normalize Fingerprint Distance 

structs = np.vstack([average_soap_LDA, average_soap_HDA, average_soap_LDA,average_soap_Isd])
distance = squareform(pdist(structs))
# Calculate the minimum and maximum distances
min_distance = np.min(distance)
max_distance = np.max(distance)

# Normalize the distance matrix
normalized_fingerprint_distance = 1- ((distance - min_distance) / (max_distance - min_distance))

# Spider Plot for Fingerprint Visualization 

# Create a list to hold the data for all three structures
data = {
    "Structure1": ["Starting Structure", "Starting Structure", "Starting Structure"],
    "Structure2": ["HDA", "LDA", "cu_Si"],
    "Normalized_Fingerpring_Distance": [fingerprint_distance_HDA, fingerprint_distance_LDA, fingerprint_distance_cuSi]
}

# Create the circle plot
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 6))
# Set the background color of the spider plot
#fig.patch.set_facecolor('lightgray')  # Set the desired background color
ax.set_facecolor('whitesmoke')
ax.xaxis.grid(True, linestyle='-', color='black')
ax.yaxis.grid(True, linestyle='--', color='black') 

# Number of vertices for the circle plot
num_vertices = len(df_start)

# Create a list of angles for the vertices
angles = [n / float(num_vertices) * 2 * np.pi for n in range(num_vertices)]
angles += angles[:1]  # To close the shape

# Plot the vertex values for the starting structure
values_start = df_start["Normalized_Fingerpring_Distance"].tolist()
values_start += values_start[:1]  # To close the shape
ax.plot(angles, values_start, linewidth=2, linestyle="solid", marker="o", label="Start", color="blue")
ax.fill(angles, values_start, alpha=0.25)
# Set the labels for the vertices
ax.set_xticks(angles[:-1])
ax.set_xticklabels(df_start["Structure2"].tolist(),fontsize = 12)


# Move the legend outside the plot
#ax.legend(loc='upper left',prop={'size': 13.5}, bbox_to_anchor=(0.9, 1.0))
$plt.savefig('{}/sl_SOAP_figs/1000_sl.png'.format(os.getcwd()), format='png', dpi=300, bbox_inches='tight')
# Display the plot
plt.show()


