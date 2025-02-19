# my_3d_descriptors/descriptors/approach_a.py

from .base_descriptor import Base3DDescriptor
from ..utils.chem_utils import _compute_MC2, _compute_fraction_sp3, _compute_fraction_carbon

class ApproachA(Base3DDescriptor):
    """ This is the first approach for the 3D descriptor. 
    Approach A:
    The 3D coordinates will be the following: 
      x = MC2 (custom complexity score). For more information on MC2 see compute_MC2() method. 
      y = fraction of sp3-hybridized atoms
      z = fraction of carbon atoms

    Args:
        Base3DDescriptor (_type_): _description_
    """

    def compute(self, mol) -> tuple[float, float, float]:
        """_summary_

        Args:
            mol (_type_): _description_

        Returns:
            _type_: _description_
        """
        x = _compute_MC2(mol)          # returns float
        y = _compute_fraction_sp3(mol) # returns float
        z = _compute_fraction_carbon(mol)  # returns float

        return (x, y, z)
