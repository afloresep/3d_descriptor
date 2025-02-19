# my_3d_descriptors/descriptors/approach_b.py

from math import log
from .base_descriptor import Base3DDescriptor
from rdkit.Chem import Descriptors
from utils.chem_utils import _get_heavy_atom_count, _compute_fraction_sp3

class ApproachB(Base3DDescriptor):
    """
    Approach B:
      x = log(HAC)
      y = TPSA / HAC
      z = fraction of sp3-hybridized atoms
    """

    def compute(self, mol):
        hac = _get_heavy_atom_count(mol)
        if hac == 0:
            return (0.0, 0.0, 0.0)  # or handle edge cases differently
        
        x = log(hac)
        tpsa = Descriptors.TPSA(mol)
        y = tpsa / hac
        z = _compute_fraction_sp3(mol)

        return (x, y, z)