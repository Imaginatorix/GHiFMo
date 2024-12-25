# Import necessary libraries
from rdkit import Chem
from rdkit.Chem import AllChem
import psi4

# Step 1: Define the drug-like molecule as a SMILES string
smiles = "CCOc1ccc2nc(S(N)(=O)=O)sc2c1"  # Example: Sulfathiazole, a drug-like molecule

# Step 2: Convert SMILES to 3D geometry using RDKit
def smiles_to_geometry(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # Add hydrogens
    AllChem.EmbedMolecule(mol)  # Generate initial 3D coordinates
    AllChem.UFFOptimizeMolecule(mol)  # Optimize geometry with UFF force field

    # Extract atomic symbols and 3D coordinates
    conf = mol.GetConformer()
    geometry = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        geometry.append(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}")
    return "\n".join(geometry)

geometry = smiles_to_geometry(smiles)

# Step 3: Define and run a DFT calculation with Psi4
def calculate_homo_lumo(geometry):
    # Initialize Psi4 molecule
    psi4.geometry(f"""
    {geometry}
    """)

    # Set Psi4 options
    psi4.set_options({
        "basis": "sto-3g",  # Basis set
        "scf_type": "df",  # Density fitting for speed
        "reference": "rhf",  # Restricted HF for closed-shell systems
    })

    # Perform DFT calculation with B3LYP
    energy, wavefunction = psi4.energy("b3lyp", return_wfn=True)

    # Extract orbital energies
    orbital_energies = wavefunction.epsilon_a().to_array()
    num_electrons = wavefunction.nalpha() + wavefunction.nbeta()
    homo_index = num_electrons // 2 - 1  # Index of HOMO
    lumo_index = homo_index + 1  # Index of LUMO

    homo_energy = orbital_energies[homo_index] * 27.2114  # Convert Hartree to eV
    lumo_energy = orbital_energies[lumo_index] * 27.2114  # Convert Hartree to eV
    homo_lumo_gap = lumo_energy - homo_energy

    return homo_energy, lumo_energy, homo_lumo_gap

homo_energy, lumo_energy, homo_lumo_gap = calculate_homo_lumo(geometry)

# Step 4: Print the results
print(f"HOMO energy: {homo_energy:.4f} eV")
print(f"LUMO energy: {lumo_energy:.4f} eV")
print(f"HOMO-LUMO gap: {homo_lumo_gap:.4f} eV")
