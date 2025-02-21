# API -class based- to compute the fingerprints
import logging
from typing import Optional, List, Union
import os
from multiprocessing import Pool
from typing import List, Tuple
import pandas as pd
from rdkit import Chem
from tqdm import tqdm
from descriptors.a_method import ApproachA  


def _compute_single_smiles(smiles: str) -> Tuple[float, float, float]:
    """
    Helper function for multiprocessing. Takes smiles and returns the fingperprint (tuple)

    Args:
        smiles (str): The SMILES string representing a molecule.

    Returns:
        tuple[float, float, float]: A tuple with (x, y, z) 3D descriptor values.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (0.0, 0.0, 0.0)

    threefp = ApproachA()

    return threefp.compute(mol)


class FingerprintCalculator:
    """
    A class that provides methods to compute fingerprints (3D descriptors) for a
    set of SMILES strings, both in-memory and in batches from a file.
    """

    def compute_fingerprints(self, smiles_list: List[str], num_processes: int = os.cpu_count()) -> List[Tuple[float, float, float]]:
        """Compute 3D descriptor values for a list of SMILES strings.

        Args:
            smiles_list (List[str]): A list of SMILES strings to compute descriptors for.
            num_processes (int): Number of parallel processes to use. Default is all cpu

        Returns:
            List[tuple[float, float, float]]: List of (x, y, z) descriptor tuples.
        """
        with Pool(processes=num_processes) as pool:
            results = pool.map(_compute_single_smiles, smiles_list)
        return results


    def calculate_fingerprints(
        self,
        smiles_list: Optional[List[str]] = None,
        input_file: Optional[str] = None,
        batch_size: Optional[int] = None,
        output_file: Optional[str] = None,
        num_processes: int = os.cpu_count(), 
    ) -> Union[List[Tuple[float, float, float]], None]:
        """
        Depending on the arguments, compute fingerprints either:
          1) In memory for a given list of SMILES (if `smiles_list` is provided).
          2) For an entire file at once (if `input_file` is provided and `batch_size` is None).
          3) In batches from a file (if `input_file` is provided and `batch_size` > 0).

        Args:
            smiles_list (Optional[List[str]]): A list of SMILES strings in memory.
            input_file (Optional[str]): Path to input file containing SMILES.
            batch_size (Optional[int]): Number of SMILES to process in each batch.
                                        If None or 0, process the entire file at once.
            output_file (Optional[str]): Output file (if None, returns result in code).
            output_format (str): One of {"csv", "parquet"}, specifying output format.
                                 If no `output_file` is given, this is ignored (return data in code).
            num_processes (int): Number of parallel processes to use.

        Returns:
            Union[List[Tuple[float, float, float]], None]:
                - A list of fingerprint tuples if no output file is given.
                - None if an output file is written (since data is saved to disk).
        """
        if output_file.lower().endswith('.csv'):
            output_format = 'csv'
        elif output_file.lower().endswith('.parquet'):
            output_format = 'parquet'
        else:
            raise ValueError(f"Output file format should be 'parquet, csv', {output_file.lower()}")

        # If the user provided an in-memory list (no file):
        if smiles_list is not None and input_file is None:
            results = self.compute_fingerprints(smiles_list, num_processes=num_processes)
            if output_file is None:
                return results 
            else:
                # Convert results to a DataFrame to save
                x, y, z = zip(*results)
                df = pd.DataFrame({
                    "smiles":smiles_list,
                    "x": x,
                    "y": y,
                    "z": z,
                })
                self._save_results(df, output_file)
                return None

        #  If the user provided an input file but NO batch size
        #    or a batch_size = 0 or None => load entire file at once.
        if input_file is not None and (not batch_size or batch_size <= 0):
            # Load entire file (assuming one SMILES per line here)
            with open(input_file, "r") as f:
                smiles_list = [line.strip() for line in f if line.strip()]

            results = self.compute_fingerprints(smiles_list, num_processes=num_processes)
            if output_file is None:
                return results
            else:
                x, y, z = zip(*results)
                df = pd.DataFrame({
                    "smiles":smiles_list,
                    "x": x,
                    "y": y,
                    "z": z,
                })
                self._save_results(df, output_file)
                return None

        # If the user provided an input file AND a positive batch size => chunk process
        if input_file and batch_size and batch_size > 0:
            # process in chunks and save each chunk to disk with an index
            chunk_counter = 1
            partial_smiles = []
            file_basename, file_ext = os.path.splitext(os.path.basename(input_file))

            with open(input_file, "r") as f:
                total_lines = sum(1 for line in f)

            with open(input_file, "r") as f, tqdm(total=total_lines, desc="processing smiles", unit="mol") as pbar:
                for line in f:
                    smiles = line.strip()
                    if not smiles:
                        continue
                    partial_smiles.append(smiles)

                    if len(partial_smiles) == batch_size:
                        results = self.compute_fingerprints(partial_smiles, num_processes=num_processes)
                        #unpack 
                        x, y, z = zip(*results)
                        df = pd.DataFrame({
                            "smiles":partial_smiles,
                            "x": x,
                            "y": y,
                            "z": z,
                        })
                        chunk_file_name = f"{file_basename}-fp_{chunk_counter}{'.csv' if output_format=='csv' else '.parquet'}"
                        self._save_results(df, chunk_file_name, output_format)
                        partial_smiles.clear()
                        chunk_counter += 1
                        pbar.update(batch_size)

            # Process any remaining SMILES
            if partial_smiles:
                results = self.compute_fingerprints(partial_smiles, num_processes=num_processes)
                x, y, z = zip(*results)
                df = pd.DataFrame({
                    "smiles":partial_smiles,
                    "x": x,
                    "y": y,
                    "z": z,
                })
                chunk_file_name = f"{file_basename}-fp_{chunk_counter}{'.csv' if output_format=='csv' else '.parquet'}"
                self._save_results(df, chunk_file_name, output_format)

            return None

        # Fallback if no condition was met (should rarely happen):
        raise ValueError("Invalid argument combination. Provide either `smiles_list` or `input_file`.")

    def _save_results(self, df: pd.DataFrame, output_file: str) -> None:
        """Internal helper method for saving data in the desired format."""
        # Check the output name to decide the output format 
        if output_file.lower().endswith("csv"):
            df.to_csv(output_file, index=False)
        if output_file.lower().endswith("parquet"):
            df.to_parquet(output_file, index=False)