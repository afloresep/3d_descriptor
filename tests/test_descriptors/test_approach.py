
from rdkit import Chem

#######
# Test to see which method returns the best count of carboxyl like group method
#######

# Method 1: Simple SMARTS pattern
# -----------------------------
def count_carboxyl_like_groups_method1(mol: Chem.Mol) -> int:
    """
    Count groups using the simple SMARTS pattern: [#7,#8]C(=O)
    """
    pattern = Chem.MolFromSmarts("[#7,#8]C(=O)")
    matches = mol.GetSubstructMatches(pattern)
    return len(matches)

# Method 2: Specific carbonyl + neighbor check (using [CX3]=[OX1])
CARBONYL_PATTERN = Chem.MolFromSmarts("[CX3]=[OX1]")

def count_carboxyl_like_groups_method2(mol: Chem.Mol) -> int:
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

# Method 3: Broader SMARTS for C=O and neighbor check (using [C,c]=O)
def count_carboxyl_like_groups_method3(mol: Chem.Mol) -> int:
    """
    Count carboxyl-like groups using a broader SMARTS for the carbonyl and then 
    checking for an attached N or O (besides the carbonyl oxygen).
    """
    if mol is None:
        raise ValueError("Invalid SMILES string.")

    # Define the SMARTS pattern for the carbonyl group (C=O)
    pattern_C_O = Chem.MolFromSmarts('[C,c]=O')

    # Find all carbonyl groups in the molecule
    matches_C_O = mol.GetSubstructMatches(pattern_C_O)

    # Initialize the count of C=O with either N or O (excluding the bonded oxygen) as a neighbor
    count_with_N_or_O = 0

    for match in matches_C_O:
        carbon_idx, oxygen_idx = match  # Get the carbon and oxygen atom indices in the C=O group

        # Get neighboring atoms of the carbon atom
        neighbors = mol.GetAtomWithIdx(carbon_idx).GetNeighbors()

        # Exclude the oxygen atom from the C=O group when checking neighbors
        valid_neighbors = [neighbor for neighbor in neighbors if neighbor.GetIdx() != oxygen_idx]

        # Check if any remaining neighbor is either nitrogen (atomic number 7) or oxygen (atomic number 8)
        if any(neighbor.GetAtomicNum() in [7, 8] for neighbor in valid_neighbors):
            count_with_N_or_O += 1

    return count_with_N_or_O

def test_carboxyl_like_group_methods():
    test_smiles = {
        "Acetic Acid": "CC(=O)O",
        "Ethyl Acetate": "CC(=O)OCC",
        "Acetamide": "CC(=O)N",
        "Urea": "NC(=O)N",
        "Carbamate": "COC(=O)N",
        "Carbonate": "OC(=O)OC",
        "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "Propionic Acid": "CCC(=O)O",
        "Butyric Acid": "CCCC(=O)O",
        "Benzoic Acid": "C1=CC=C(C=C1)C(=O)O",
        "Formic Acid": "O=C(O)",
        "Succinic Acid": "O=C(O)CC(=O)O",
        "Maleic Acid": "O=C(O)/C=C\\C(=O)O",
        "Fumaric Acid": "O=C(O)/C=C/C(=O)O",
        "Oxalic Acid": "O=C(O)C(=O)O",
        "Pyruvic Acid": "CC(=O)C(=O)O",
        "Lactic Acid": "CC(O)C(=O)O",
        "Citric Acid": "C(C(=O)O)(CC(=O)O)C(O)(C(=O)O)",
        "Salicylic Acid": "OC1=CC=CC=C1C(=O)O",
        "Benzamide": "C1=CC=C(C=C1)C(=O)N",
        "Acetanilide": "CC(=O)NC1=CC=CC=C1",
        "Propionamide": "CCC(=O)N",
        "Butyramide": "CCCC(=O)N",
        "Nicotinamide": "NC(=O)c1cccnc1",
        "Formamide": "C(=O)N",
        "N-Methylformamide": "C(=O)NC",
        "N,N-Dimethylformamide": "C(=O)N(C)C",
        "Dimethylacetamide": "CC(=O)N(C)C",
        "N-Methylacetamide": "CC(=O)NC",
        "Methyl Benzoate": "C1=CC=C(C=C1)C(=O)OC",
        "Ethyl Benzoate": "C1=CC=C(C=C1)C(=O)OCC",
        "Isopropyl Acetate": "CC(=O)OC(C)C",
        "Butyl Acetate": "CC(=O)OCCCC",
        "Methyl Formate": "COC(=O)",
        "Ethyl Formate": "C(=O)OCC",
        "Propyl Acetate": "CC(=O)OCCC",
        "Benzyl Acetate": "CC(=O)OCC1=CC=CC=C1",
        "Succinimide": "O=C1NC(=O)C1",
        "Phthalimide": "O=C1c2ccccc2C(=O)N1",
        "Maleimide": "O=C1C=CN1",
        "Acetic Anhydride": "CC(=O)OC(=O)C",
        "Succinic Anhydride": "O=C1OC(=O)CC1",
        "Maleic Anhydride": "O=C1OC(=O)C=C1",
        "Phthalic Anhydride": "O=C1OC(=O)c2ccccc2C1",
        "Glutamic Acid": "NC(CCC(=O)O)C(=O)O",
        "Aspartic Acid": "NC(CC(=O)O)C(=O)O",
        "Glycine": "NCC(=O)O",
        "Alanine": "CC(N)C(=O)O",
        "Valine": "CC(C)C(N)C(=O)O",
        "p-Hydroxybenzoic Acid": "C1=CC(=CC=C1O)C(=O)O",
        "CHEBI:59968": "OCc1ccc(C=O)n1CCCCCC(O)=O", 
        "CHEBI:140643": "C=1(C=CC(=CC1)[N+]([O-])=O)COP(CCCCCC(O)=O)(=O)O", 
        "CHEBI:192408": "BrCC([O-])=O", 
        "CHEBI:59213 ": "[H][C@]12SCC(CSc3nnnn3C)=C(N1C(=O)[C@@]2([H])NC(=O)[C@H](NC(=O)c1cnc(C)cc1O)c1ccc(O)cc1)C(O)=O",
        "CHEBI:210250": "O=C(O)CC1=C(CC[C@H]2[C@]1(CCC[C@]2(CO)C)C)CC(=CCO)C",
        "CHEBI:80036": "CC[C@H](C)[C@H](NC(=O)[C@H]1CCCCN1C)C(=O)N(COC(=O)CC(C)C)[C@H](C[C@@H](OC(C)=O)c1nc(cs1)C(=O)N[C@H](C[C@H](C)C(O)=O)Cc1ccccc1)C(C)C"

    }
    print("Comparing carboxyl-like group counts:")
    print("{:<20} {:>8} {:>8} {:>8}".format("Molecule", "Method1", "Method2", "Method3"))
    print("-"*50)
    
    for name, smiles in test_smiles.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Could not parse SMILES for {name}")
            continue
        count1 = count_carboxyl_like_groups_method1(mol)
        count2 = count_carboxyl_like_groups_method2(mol)
        count3 = count_carboxyl_like_groups_method3(mol)
        print("{:<20} {:>8} {:>8} {:>8}".format(name, count1, count2, count3))

if __name__ == "__main__":
    test_carboxyl_like_group_methods()
    """
    Should return something like
    Molecule              Method1  Method2  Method3
    --------------------------------------------------
    Acetic Acid                 1        1        1
    Ethyl Acetate               1        1        1
    Acetamide                   1        1        1
    Urea                        2        1        1
    Carbamate                   2        1        1
    Carbonate                   2        1        1
    Aspirin                     2        2        2
    ...
    """
