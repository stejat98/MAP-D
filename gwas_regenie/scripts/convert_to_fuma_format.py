#!/usr/bin/env python3
"""
Convert REGENIE summary statistics files to FUMA-compatible format.

FUMA Requirements:
- Required columns: CHR, BP, SNP, A1 (effect allele), A2 (non-effect allele), P, BETA, SE, N
- Optional: A1FREQ, OR
- Files must be gzipped if >600MB
- Space or tab separated
"""

import gzip
import csv
import sys
import argparse
from pathlib import Path
import math


def convert_to_fuma_format_row(row, header_map):
    """
    Convert a single row from REGENIE format to FUMA format
    
    Parameters:
    -----------
    row : dict
        Dictionary with REGENIE column names as keys
    header_map : dict
        Mapping from REGENIE column names to FUMA column names
    """
    fuma_row = {}
    
    # Map required columns
    fuma_row['CHR'] = row.get('CHROM', row.get('CHR'))
    fuma_row['BP'] = row.get('GENPOS', row.get('BP'))
    fuma_row['SNP'] = row.get('ID', row.get('SNP'))
    fuma_row['A1'] = row.get('ALLELE1', row.get('A1'))  # Effect allele
    fuma_row['A2'] = row.get('ALLELE0', row.get('A2'))  # Non-effect allele
    
    # Calculate P from LOG10P if needed
    if 'P' in row:
        p_val = float(row['P'])
    elif 'LOG10P' in row:
        log10p = float(row['LOG10P'])
        p_val = 10 ** (-log10p)
    else:
        raise ValueError("Missing P or LOG10P column")
    
    fuma_row['P'] = p_val
    fuma_row['BETA'] = row.get('BETA')
    fuma_row['SE'] = row.get('SE')
    fuma_row['N'] = row.get('N')
    
    # Optional columns
    if 'A1FREQ' in row:
        fuma_row['A1FREQ'] = row['A1FREQ']
    if 'OR' in row:
        fuma_row['OR'] = row['OR']
    
    return fuma_row


def convert_file(input_path, output_path=None, check_size=True):
    """
    Convert a single REGENIE file to FUMA format
    
    Parameters:
    -----------
    input_path : str
        Path to input .regenie.gz file
    output_path : str
        Path to output file (default: input_path with .fuma.gz extension)
    check_size : bool
        Check if file size is acceptable for FUMA (<600MB uncompressed)
    """
    input_path = Path(input_path)
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    if output_path is None:
        # Create output filename: replace .regenie.gz with .fuma.gz
        output_path = input_path.parent / input_path.name.replace('.regenie.gz', '.fuma.gz')
    else:
        output_path = Path(output_path)
    
    print(f"\n{'='*60}")
    print(f"Converting: {input_path.name}")
    print(f"{'='*60}")
    
    # Read the REGENIE file
    print(f"Reading file...")
    
    # Open input file
    with gzip.open(input_path, 'rt') as f_in:
        # Read header
        header_line = f_in.readline().strip()
        regenie_headers = header_line.split()
        
        # Verify required columns exist
        required_regenie_cols = ['CHROM', 'GENPOS', 'ID', 'ALLELE0', 'ALLELE1', 'BETA', 'SE', 'N']
        missing_cols = [col for col in required_regenie_cols if col not in regenie_headers]
        
        if missing_cols:
            raise ValueError(f"Missing required columns in input file: {missing_cols}")
        
        if 'LOG10P' not in regenie_headers and 'P' not in regenie_headers:
            raise ValueError("Missing LOG10P or P column in input file")
        
        # Create header map
        header_map = {h: h for h in regenie_headers}
        
        # Define FUMA output columns
        fuma_headers = ['CHR', 'BP', 'SNP', 'A1', 'A2', 'P', 'BETA', 'SE', 'N']
        if 'A1FREQ' in regenie_headers:
            fuma_headers.append('A1FREQ')
        if 'OR' in regenie_headers:
            fuma_headers.append('OR')
        
        # Read and convert rows
        print("Converting rows...")
        rows_processed = 0
        rows_written = 0
        
        # Open output file
        with gzip.open(output_path, 'wt') as f_out:
            writer = csv.DictWriter(f_out, fieldnames=fuma_headers, delimiter='\t')
            writer.writeheader()
            
            # Process each row
            reader = csv.DictReader(f_in, fieldnames=regenie_headers, delimiter=' ')
            
            for row in reader:
                rows_processed += 1
                
                # Skip empty rows
                if not row.get('CHROM') or not row.get('ID'):
                    continue
                
                try:
                    # Convert row
                    fuma_row = convert_to_fuma_format_row(row, header_map)
                    
                    # Validate values
                    p_val = float(fuma_row['P'])
                    n_val = float(fuma_row['N'])
                    se_val = float(fuma_row['SE'])
                    
                    # Filter invalid rows
                    if p_val <= 0 or p_val > 1 or n_val <= 0 or se_val <= 0:
                        continue
                    
                    # Write row
                    writer.writerow(fuma_row)
                    rows_written += 1
                    
                    # Progress update every 1M rows
                    if rows_processed % 1000000 == 0:
                        print(f"  Processed {rows_processed:,} rows, written {rows_written:,} valid rows...")
                
                except (ValueError, KeyError) as e:
                    # Skip rows with conversion errors
                    continue
        
        print(f"Loaded {rows_processed:,} variants")
        print(f"After conversion: {rows_written:,} variants")
    
    # Check file size
    if check_size:
        final_size_mb = output_path.stat().st_size / (1024**2)
        print(f"Output file size: {final_size_mb:.1f} MB (compressed)")
        
        if final_size_mb > 600:
            print(f"⚠️  WARNING: Compressed file size exceeds 600MB (FUMA limit)")
            print(f"   You may need to filter variants (e.g., by MAF or P-value)")
    
    print(f"✓ Conversion complete!")
    print(f"  Output: {output_path}")
    print(f"  Columns: {', '.join(fuma_headers)}")
    
    return output_path


def main():
    parser = argparse.ArgumentParser(
        description='Convert REGENIE summary statistics to FUMA-compatible format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert a single file
  python convert_to_fuma_format.py BMI_diabetes_all_chr.regenie.gz
  
  # Convert all files in results directory
  python convert_to_fuma_format.py --all
  
  # Convert with custom output name
  python convert_to_fuma_format.py input.regenie.gz -o output.fuma.gz
        """
    )
    
    parser.add_argument('input_file', nargs='?', help='Path to input .regenie.gz file')
    parser.add_argument('--all', action='store_true',
                       help='Convert all .regenie.gz files in results/ directory')
    parser.add_argument('--output', '-o', help='Output file path (optional)')
    parser.add_argument('--results-dir', default='results',
                       help='Results directory (default: results)')
    parser.add_argument('--no-size-check', action='store_true',
                       help='Skip file size checking')
    
    args = parser.parse_args()
    
    if args.all:
        # Convert all files in results directory
        results_dir = Path(args.results_dir)
        if not results_dir.exists():
            print(f"Error: Results directory not found: {results_dir}")
            sys.exit(1)
        
        regenie_files = list(results_dir.glob('*.regenie.gz'))
        
        if not regenie_files:
            print(f"No .regenie.gz files found in {results_dir}")
            sys.exit(1)
        
        print(f"Found {len(regenie_files)} files to convert:")
        for f in regenie_files:
            print(f"  - {f.name}")
        
        # Convert each file
        converted_files = []
        for input_file in regenie_files:
            try:
                output_file = convert_file(input_file, check_size=not args.no_size_check)
                converted_files.append(output_file)
            except Exception as e:
                print(f"✗ Error converting {input_file.name}: {e}")
                import traceback
                traceback.print_exc()
                continue
        
        print(f"\n{'='*60}")
        print(f"Conversion Summary")
        print(f"{'='*60}")
        print(f"Successfully converted {len(converted_files)} files:")
        for f in converted_files:
            print(f"  ✓ {f.name}")
    
    elif args.input_file:
        # Convert single file
        try:
            convert_file(args.input_file, args.output, check_size=not args.no_size_check)
        except Exception as e:
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
