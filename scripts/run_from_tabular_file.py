from tqdm import tqdm
import os
import argparse
from fingerprints.fingerprint_calculator import FingerprintCalculator
from fingerprints.fingerprint_calculator import FingerprintCalculator


##############################################################################
# Script to calculate and save fingerprints from a file with tabular (\t) data
##############################################################################

def process_single_file():
    """This script will run a single file. Progress bar will be updated
    by number of molecules processed. 
    """
    parser = argparse.ArgumentParser(
        description="Compute 3D descriptors for SMILES data."
    )
    parser.add_argument(
        "--input-file", 
        type=str, 
        help="Path to input file (one SMILES per line)."
    )
    parser.add_argument(
        "--batch-size", 
        type=int, 
        default=0,
        help="If > 0, process in chunks of this size. If 0, load entire file at once."
    )
    parser.add_argument(
        "--num-processes",
        type=int,
        default=os.cpu_count(),
        help="Number of parallel processes for fingerprint calculation."
    )
    parser.add_argument(
        "--smiles-index", 
        type=int,
        default=0,
        help="Position where the smiles is on the tabular input file. Default is 0"
    )

    args = parser.parse_args()

    calculator = FingerprintCalculator()

    # Define some helpful variables
    batch_count = 0
    smiles_list =[]

    # Get total number of lines
    # This is slow for big files 
    with open(args.input_file, 'r') as file:
        total_lines = sum(1 for _ in file)

    # Iterative process for the file
    with open(args.input_file, 'r') as file, tqdm(desc="Processing file", unit="mol", total=total_lines) as pbar:
        for line in file:
            smiles= line.split("\t")[args.smiles_index] # this should specify the position of the smiles strings in your tabular data -usally 0-. 
            smiles_list.append(smiles)
            if len(smiles_list) == args.batch_size:
                calculator.calculate_fingerprints(smiles_list=smiles_list, output_file=f"/mnt/10tb_hdd/new_3fp/output_{batch_count}.parquet")
                batch_count +=1
                smiles_list =[]
                pbar.update(args.batch_size)

def process_multiple_file():
    """This script will run a single file. Progress bar will be updated
    by number of molecules processed. 
    """
    parser = argparse.ArgumentParser(
        description="Compute 3D descriptors for SMILES data."
    )
    parser.add_argument(
        "--input-folder", 
        type=str, 
        help="Path to input folder, each file in there must have one SMILES per line."
    )
    parser.add_argument(
        "--batch-size", 
        type=int, 
        default=0,
        help="If > 0, process in chunks of this size. If 0, load entire file at once."
    )
    parser.add_argument(
        "--num-processes",
        type=int,
        default=os.cpu_count(),
        help="Number of parallel processes for fingerprint calculation."
    )
    parser.add_argument(
        "--smiles-index", 
        type=int,
        default=0,
        help="Position where the smiles is on the tabular input file. Default is 0"
    )

    args = parser.parse_args()

    calculator = FingerprintCalculator()

    # Define some helpful variables
    batch_count = 0
    smiles_list =[]

    with tqdm(desc="Processing files", unit="file", total=len(os.listdir(args.input_folder))) as pbar:
        for input_file in os.listdir(args.input_folder):         # Iterative process for the file
            with open(os.path.join(os.path.abspath(args.input_folder), input_file ), 'r') as file:
                for line in file:
                    # smiles= line.split("\t")[args.smiles_index] # this should specify the position of the smiles strings in your tabular data -usally 0-. 
                    smiles_list.append(line)
                    if len(smiles_list) == args.batch_size:
                        calculator.calculate_fingerprints(smiles_list=smiles_list, output_file=f"{input_file.split('.')[0]}_output_{batch_count}.parquet")
                        batch_count +=1
                        smiles_list =[]
                if len(smiles_list) != 0:
                # Save remaining lines that are less than batch size 
                    calculator.calculate_fingerprints(smiles_list=smiles_list, output_file=f"{input_file.split('.')[0]}_output_{batch_count}.parquet")
            pbar.update()

if __name__=="__main__":
    process_single_file()