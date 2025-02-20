from rdkit import Chem

def _compute_MC2(mol):
    """Computes the MC2 descriptor value. 
    MC2= HAC - DV - 2x(X-C=O),  
    Args:
        mol (Chem.Mol): RDKit molecule object.
    
    Returns:
        int: MC2 complexity score.
    Where 
		- HAC = heavy atom count,
        - DV = number of divalent nodes (like bridging atoms in the molecular graph),
		- (X-C=O) = number of carbonyl + heteroatom “carboxyl-like” groups (amides, esters, carboxylic acids, etc.).
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

def _compute_MC2(mol: Chem.Mol) -> int:
    """Computes the MC2 descriptor value.
    
    MC2 = HAC - DV - 2x(X-C=O)
  
    """
    hac = _get_heavy_atom_count(mol)
    dv = _count_divalent_nodes(mol)
    x_c_o = _count_carboxyl_like_groups(mol)
    mc2 = hac - dv - 2 * x_c_o
    return mc2




def _get_heavy_atom_count(mol: Chem.Mol) -> int:
    """Computes the heavy atom count (HAC).
    
    Args:
        mol (Chem.Mol): RDKit molecule object.
    
    Returns:
        int: Number of heavy atoms in the molecule.
    """
    return mol.GetNumHeavyAtoms()

def _compute_fraction_sp3(mol: Chem.Mol) -> float:
    """Computes the fraction of sp3-hybridized atoms.
    
    Args:
        mol (Chem.Mol): RDKit molecule object.
    
    Returns:
        float: Fraction of sp3-hybridized atoms among heavy atoms.
    """
    sp3_count = sum(1 for atom in mol.GetAtoms() if atom.GetHybridization() == Chem.HybridizationType.SP3)
    hac = _get_heavy_atom_count(mol)
    return sp3_count / hac if hac > 0 else 0.0


def _compute_fraction_carbon(mol: Chem.Mol) -> float:
    """Computes the fraction of carbon atoms in the molecule.
    
    Args:
        mol (Chem.Mol): RDKit molecule object.
    
    Returns:
        float: Fraction of carbon atoms among heavy atoms.
    """
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    hac = _get_heavy_atom_count(mol)
    return carbon_count / hac if hac > 0 else 0.0

def _count_divalent_nodes(mol: Chem.Mol) -> int:
    """Counts the number of divalent nodes in the molecule.
    
    Args:
        mol (Chem.Mol): RDKit molecule object.
    
    Returns:
        int: Number of atoms with exactly two neighbors.

    In a more easy to read format: 
    divalent_nodes = 0
    for atom in mol.GetAtoms():
        if atom.GetDegree() == 2:
            divalent_nodes += 1 
    return divalent_nodes
    """
    return sum(1 for atom in mol.GetAtoms() if atom.GetDegree() == 2) 

def _count_carboxyl_like_groups(mol: Chem.Mol) -> int:
    """Counts the number of carboxyl-like groups in the molecule.
    
    A carboxyl-like group is defined as a carbonyl carbon attached to at least one nitrogen or oxygen atom.
    
    Args:
        mol (Chem.Mol): RDKit molecule object.
    
    Returns:
        int: Number of carboxyl-like groups.
    """
    # This should be somehow chagned to not be computed every time the function is called. 
    CARBONYL_PATTERN = Chem.MolFromSmarts("[CX3]=[OX1]")

    matches_c_o = mol.GetSubstructMatches(CARBONYL_PATTERN)
    count_with_N_or_O = 0
    counted_carbons = set()
    
    for match in matches_c_o:
        if len(match) != 2:
            continue
        carbon_idx, oxygen_idx = match
        if carbon_idx in counted_carbons:
            continue
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        neighbors = carbon_atom.GetNeighbors()
        valid_neighbors = [nbr for nbr in neighbors if nbr.GetIdx() != oxygen_idx]
        if any(nbr.GetAtomicNum() in [7, 8] for nbr in valid_neighbors):
            count_with_N_or_O += 1
            counted_carbons.add(carbon_idx)
    
    return count_with_N_or_O