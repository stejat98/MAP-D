#!/usr/bin/env python3
"""
Script to read REGENIE summary statistics files
Compatible with UK Biobank GWAS format standards
"""

import gzip
import pandas as pd
import sys
import argparse
from pathlib import Path


def read_regenie_summary_stats(file_path, add_pvalue=True, convert_to_ukbb_format=False):
    """
    Read REGENIE summary statistics file
    
    Parameters:
    -----------
    file_path : str
        Path to the .regenie.gz file
    add_pvalue : bool
        Convert LOG10P to P-value column
    convert_to_ukbb_format : bool
        Convert column names to UK Biobank format
    
    Returns:
    --------
    pd.DataFrame
        Summary statistics dataframe
    """
    print(f"Reading file: {file_path}")
    
    # Read gzipped file
    df = pd.read_csv(file_path, sep=' ', compression='gzip', low_memory=False)
    
    print(f"Loaded {len(df):,} variants")
    
    # Convert LOG10P to P-value if requested
    if add_pvalue and 'LOG10P' in df.columns:
        df['P'] = 10 ** (-df['LOG10P'])
        print("Added P-value column (from LOG10P)")
    
    # Convert to UK Biobank format if requested
    if convert_to_ukbb_format:
        df = convert_to_ukbb_format(df)
    
    return df


def convert_to_ukbb_format(df):
    """
    Convert REGENIE format to UK Biobank format
    
    UK Biobank standard columns:
    CHR, BP, SNP, A1, A2, FRQ, N, BETA, SE, P
    """
    df_ukbb = df.rename(columns={
        'CHROM': 'CHR',
        'GENPOS': 'BP',
        'ID': 'SNP',
        'ALLELE1': 'A1',  # Effect allele
        'ALLELE0': 'A2',  # Reference allele
        'A1FREQ': 'FRQ'
    })
    
    # Select standard columns (P should already exist if add_pvalue=True)
    standard_cols = ['CHR', 'BP', 'SNP', 'A1', 'A2', 'FRQ', 'N', 'BETA', 'SE', 'P']
    if 'CHISQ' in df_ukbb.columns:
        standard_cols.append('CHISQ')
    if 'LOG10P' in df_ukbb.columns:
        standard_cols.append('LOG10P')
    
    df_ukbb = df_ukbb[[col for col in standard_cols if col in df_ukbb.columns]]
    
    print("Converted to UK Biobank format")
    return df_ukbb


def filter_summary_stats(df, min_maf=0.01, max_maf=0.99, min_n=100, remove_missing=True):
    """
    Apply standard QC filters to summary statistics
    """
    n_before = len(df)
    
    if remove_missing:
        df = df.dropna(subset=['BETA', 'SE', 'P'])
    
    if 'A1FREQ' in df.columns:
        df = df[(df['A1FREQ'] >= min_maf) & (df['A1FREQ'] <= max_maf)]
    
    if 'N' in df.columns:
        df = df[df['N'] >= min_n]
    
    n_after = len(df)
    print(f"Filtered: {n_before:,} -> {n_after:,} variants (removed {n_before - n_after:,})")
    
    return df


def main():
    parser = argparse.ArgumentParser(
        description='Read REGENIE summary statistics files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python read_regenie_summary_stats.py BMI_diabetes_all_chr.regenie.gz
  python read_regenie_summary_stats.py BMI_diabetes_all_chr.regenie.gz --ukbb --filter
        """
    )
    
    parser.add_argument('file_path', help='Path to .regenie.gz file')
    parser.add_argument('--ukbb', action='store_true', 
                       help='Convert to UK Biobank format')
    parser.add_argument('--filter', action='store_true',
                       help='Apply standard QC filters')
    parser.add_argument('--output', '-o', help='Output file path (optional)')
    
    args = parser.parse_args()
    
    # Read the file
    df = read_regenie_summary_stats(
        args.file_path,
        add_pvalue=True,
        convert_to_ukbb_format=args.ukbb
    )
    
    # Apply filters if requested
    if args.filter:
        df = filter_summary_stats(df)
    
    # Print summary
    print("\n=== Summary Statistics ===")
    print(f"Total variants: {len(df):,}")
    print(f"Chromosomes: {sorted(df['CHROM'].unique())}")
    print(f"Sample size range: {df['N'].min():.0f} - {df['N'].max():.0f}")
    if 'P' in df.columns:
        print(f"P-value range: {df['P'].min():.2e} - {df['P'].max():.2e}")
    
    # Show first few rows
    print("\n=== First 5 variants ===")
    print(df.head())
    
    # Save to file if requested
    if args.output:
        df.to_csv(args.output, sep='\t', index=False, compression='gzip' if args.output.endswith('.gz') else None)
        print(f"\nSaved to: {args.output}")


if __name__ == '__main__':
    main()
