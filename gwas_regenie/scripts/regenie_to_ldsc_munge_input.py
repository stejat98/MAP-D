#!/usr/bin/env python3
"""Stream REGENIE .regenie.gz to tab-delimited SNP/A1/A2/N/P/BETA for LDSC munge_sumstats.py."""
import argparse
import gzip
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("regenie_gz", help="Input REGENIE .regenie.gz")
    p.add_argument("out_gz", help="Output .txt.gz (tab-separated with header)")
    args = p.parse_args()

    n_in = n_out = 0
    with gzip.open(args.regenie_gz, "rt") as fin, gzip.open(
        args.out_gz, "wt", compresslevel=6
    ) as fout:
        header = fin.readline()
        fout.write("\t".join(["SNP", "A1", "A2", "N", "P", "BETA"]) + "\n")
        for line in fin:
            n_in += 1
            if n_in % 2_000_000 == 0:
                print(f"  ... read {n_in:,} variant lines", file=sys.stderr)
            parts = line.split()
            if len(parts) < 12:
                continue
            try:
                log10p = float(parts[11])
            except ValueError:
                continue
            if log10p != log10p:  # NaN
                continue
            try:
                beta = float(parts[8])
                n = float(parts[6])
            except ValueError:
                continue
            pval = 10.0 ** (-log10p)
            if pval >= 1.0:
                pval = 1.0 - 1e-16
            if pval <= 0.0:
                pval = 1e-300
            snp, a1, a2 = parts[2], parts[4], parts[3]
            fout.write(
                "\t".join(
                    [
                        snp,
                        a1,
                        a2,
                        str(int(n)) if n == int(n) else str(n),
                        f"{pval:.6e}",
                        f"{beta:.8e}",
                    ]
                )
                + "\n"
            )
            n_out += 1

    print(f"Done: {n_in:,} lines read, {n_out:,} variants written", file=sys.stderr)


if __name__ == "__main__":
    main()
