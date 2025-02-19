from rdkit import Chem

def _compute_MC2(mol):
    """_summary_
    MC2= HAC - DV - 2x(X-C=O),
    where
	•	HAC = heavy atom count,
	•	DV = number of divalent nodes (like bridging atoms in the molecular graph),
	•	(X-C=O) = number of carbonyl + heteroatom “carboxyl-like” groups (amides, esters, carboxylic acids, etc.).
    	It is more nuanced than a raw size measure because it penalizes certain structural motifs 
        (e.g., bridging atoms, carbonyl-with-heteroatoms)
    """
    # Compute HAC
    hac = _get_heavy_atom_count(mol)

    # Count divalent nodes
    dv = _count_divalent_nodes(mol)

    # Count X-C=O groups (carboxyl-like)
    x_c_o = _count_carboxyl_like_groups(mol)

    # apply mc2 formula (Ye's paper)
    mc2 = hac - dv - 2*(x_c_o)

    return mc2 

def _get_heavy_atom_count(mol):
    return mol.GetNumHeavyAtoms()

def _compute_fraction_sp3(mol):
    """
    E.g. fraction of sp3-hybridized carbons among heavy atoms.
    Could be something like (sp3Carbons / HAC).
    """
    # Implementation details...
    pass

def _compute_fraction_carbon(mol):
    """
    fraction = number_of_carbon_atoms / HAC
    """
    pass

def _count_divalent_nodes(mol)-> int:
    """
    Count atoms that have exactly two neighbors.

    returns number of divalent nodes. 
    In a more easy to read format: 
    divalent_nodes = 0
    for atom in mol.GetAtoms():
        if atom.GetDegree() == 2:
            divalent_nodes += 1 
    return divalent_nodes
    """
    return sum(1 for atom in mol.GetAtoms() if atom.GetDegree() == 2) 

# Here, we use a more specific pattern for a carbonyl group.
CARBONYL_PATTERN = Chem.MolFromSmarts("[CX3]=[OX1]")

CARBONYL_PATTERN = Chem.MolFromSmarts("[CX3]=[OX1]")

def _count_carboxyl_like_groups(mol: Chem.Mol) -> int:
    """
    Count carboxyl-like groups using a specific carbonyl SMARTS and then verifying 
    that the carbonyl carbon (besides the carbonyl oxygen) is attached to at least 
    one nitrogen or oxygen.
    """
    matches_c_o = mol.GetSubstructMatches(CARBONYL_PATTERN)
    count_with_N_or_O = 0
    counted_carbons = set()  # avoid double-counting the same carbon atom

    for match in matches_c_o:
        # match is a tuple (carbon_idx, oxygen_idx)
        if len(match) != 2:
            continue
        carbon_idx, oxygen_idx = match
        if carbon_idx in counted_carbons:
            continue
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        neighbors = carbon_atom.GetNeighbors()
        # Exclude the oxygen atom in the C=O group
        valid_neighbors = [nbr for nbr in neighbors if nbr.GetIdx() != oxygen_idx]
        # Check if any remaining neighbor is N (7) or O (8)
        if any(nbr.GetAtomicNum() in [7, 8] for nbr in valid_neighbors):
            count_with_N_or_O += 1
            counted_carbons.add(carbon_idx)
    return count_with_N_or_O

