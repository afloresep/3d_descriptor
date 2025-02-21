import os
import argparse
from fingerprints.fingerprint_calculator import FingerprintCalculator


def main():
    parser = argparse.ArgumentParser(
        description="Compute 3D descriptors for SMILES data."
    )
    parser.add_argument(
        "--smiles-list", 
        nargs="+", 
        help="List of SMILES strings to compute descriptors in memory."
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
        "--output-file", 
        type=str, 
        default=None,
        help="Path to output file. If not specified, results are returned in code (printed)."
    )
    parser.add_argument(
        "--output-format",
        choices=["csv", "parquet"],
        default="csv",
        help="Output format if writing results to disk (default='csv')."
    )
    parser.add_argument(
        "--num-processes",
        type=int,
        default=os.cpu_count(),
        help="Number of parallel processes for fingerprint calculation."
    )

    args = parser.parse_args()

    calculator = FingerprintCalculator()

    if args.smiles_list and not args.input_file:
        # No file, just an in-memory list
        result = calculator.calculate_fingerprints(
            smiles_list=args.smiles_list,
            input_file=None,
            batch_size=None,  # Not applicable
            output_file=args.output_file,
            output_format=args.output_format,
            num_processes=args.num_processes
        )

        if result is not None:
            # If output_file was None, we have a list of tuples we can print or use
            for smi, desc in zip(args.smiles_list, result):
                print(f"{smi} => {desc}")

    elif args.input_file:
        # Using a file, decide how to handle batch_size
        result = calculator.calculate_fingerprints(
            smiles_list=None,
            input_file=args.input_file,
            batch_size=args.batch_size if args.batch_size > 0 else None,
            output_file=args.output_file,
            output_format=args.output_format,
            num_processes=args.num_processes
        )
        if result is not None:
            # Means we did a full file load without an output_file
            # so we have a list of descriptor tuples
            for fp_tuple in result:
                print(fp_tuple)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
