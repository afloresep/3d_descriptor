# my_3d_descriptors/utils/chem_utils.py

from rdkit import Chem

def compute_MC2(mol):
    """_summary_
    MC2= HAC - DV - 2x(X-C=O),
    where
	•	HAC = heavy atom count,
	•	DV = number of divalent nodes (like bridging atoms in the molecular graph),
	•	(X-C=O) = number of carbonyl + heteroatom “carboxyl-like” groups (amides, esters, carboxylic acids, etc.).
    	It is more nuanced than a raw size measure because it penalizes certain structural motifs 
        (e.g., bridging atoms, carbonyl-with-heteroatoms)
    """
    # 1) Compute HAC
    hac = get_heavy_atom_count(mol)

    # 2) Count divalent nodes
    dv = count_divalent_nodes(mol)

    # 3) Count X-C=O groups (carboxyl-like)
    x_c_o = count_carboxyl_like_groups(mol)

    return hac - dv - 2*x_c_o

def get_heavy_atom_count(mol):
    return mol.GetNumHeavyAtoms()

def compute_fraction_sp3(mol):
    """
    E.g. fraction of sp3-hybridized carbons among heavy atoms.
    Could be something like (sp3Carbons / HAC).
    """
    # Implementation details...
    pass

def compute_fraction_carbon(mol):
    """
    fraction = number_of_carbon_atoms / HAC
    """
    pass

def count_divalent_nodes(mol):
    """
    Count atoms that have exactly two neighbors.
    """
    pass

def count_carboxyl_like_groups(mol):
    """
    Identify patterns for carboxylic acids, esters, amides, carbamates, etc.
    Possibly use SMARTS patterns.
    """
    pass
