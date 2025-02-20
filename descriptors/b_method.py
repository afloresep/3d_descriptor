# my_3d_descriptors/descriptors/approach_b.py

from math import log
from .base_descriptor import Base3DDescriptor
from rdkit.Chem import Descriptors
from utils.chem_utils import _get_heavy_atom_count, _compute_fraction_sp3
from rdkit import Chem
class ApproachB(Base3DDescriptor):
    """Approach B:
      x = log(HAC)
      y = TPSA / HAC
      z = fraction of sp3-hybridized atoms
    
    Args:
        Base3DDescriptor (ABC): Base class for 3D descriptors.
    """

    def compute(self, mol: Chem.Mol) -> tuple[float, float, float]:
        """Computes the 3D descriptor values for the given molecule.

        Args:
            mol (Chem.Mol): RDKit molecule object.

        Returns:
            tuple[float, float, float]: A tuple containing:
                - x: log of heavy atom count (float)
                - y: TPSA divided by heavy atom count (float)
                - z: Fraction of sp3-hybridized atoms (float)
        """
        hac = _get_heavy_atom_count(mol)
        if hac == 0:
            return (0.0, 0.0, 0.0)  # or handle edge cases differently
        
        x = log(hac)
        tpsa = Descriptors.TPSA(mol)
        y = tpsa / hac
        z = _compute_fraction_sp3(mol)

        return (x, y, z)
