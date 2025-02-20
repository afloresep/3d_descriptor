from abc import ABC, abstractmethod


class Base3DDescriptor(ABC):
    """
    Abstract base class for 3D descriptors. 
    Ensures each descriptor has a compute(mol) -> (float, float, float) method. 

    Args:
        ABC (_type_): _description_
    """

    @abstractmethod
    def compute(self, mol):
        """Compute the 3D descriptor for a given moelcule

        Args:
            mol (rdkit.Chem.Mol): RDKit molecule object or a similar
            compatible object.
        Returns:
            Tuple[float, float, float]:  A tuple representing the (x, y, z) 3D
            coordinates calculated from the molecule object
        """
        pass
