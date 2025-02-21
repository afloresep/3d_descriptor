from .base_descriptor import Base3DDescriptor
from utils.chem_utils import _compute_MC2, _compute_fraction_sp3, _compute_fraction_carbon, _count_carboxyl_like_groups, _count_divalent_nodes, _get_heavy_atom_count
from rdkit import Chem
from rdkit.Chem import Descriptors
from math import log


class ApproachC(Base3DDescriptor):
    """This is the third approach for the 3D descriptor.
    Approach C: 
    The 3D coordinates will be the following: 
     x = log(MAC_2)
     y = TPSA / HAC
     z = fraction of carbon atoms

    Args:
        Base3DDescriptor (_type_): _description_

    Returns:
        _type_: _description_
    """

    def compute(self, mol:Chem.Mol) -> tuple[float, float, float]:
        hac = _get_heavy_atom_count(mol) 
        
        mac2 = _compute_MC2(mol)

        if mac2 ==0:
            mac2 = 0.00000000001 
        try: 
            x = log(mac2)
        except Exception as e:
            print(f"Raised exception {e} when doing log{mac2}")

        tpsa = Descriptors.TPSA(mol)
        y = tpsa / hac

        z = _compute_fraction_sp3(mol)

        return (x, y, z)