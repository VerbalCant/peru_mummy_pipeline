#!/usr/bin/env python3
#
# I lovehate conda.
# 
# rm -rf workflow && conda activate quast && ./binning_workflow.py   --num_samples 1000000 -o workflow SRR20458000_megahit_contigs.fasta.gz SRR20458000_spades_contigs.fasta.gz --no-checkm && conda activate checkm && ./binning_workflow.py   --num_samples 1000000 -o workflow SRR20458000_megahit_contigs.fasta.gz SRR20458000_spades_contigs.fasta.gz --no-subset --no-maxbin
# xstreme --p SRR20458000_megahit_contigs2.fasta --oc SRR20458000_megahit_xstreme --meme-p 6 --meme-searchsize 200000


import argparse
import subprocess
import os

def run_command(command):
    try:
        print(f"Running: {command}")
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {e}")
        exit(1)

def main(args):
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Subset contigs using seqtk
    megahit_subset = os.path.join(args.output_dir, 'megahit_subset')
    spades_subset = os.path.join(args.output_dir, 'spades_subset')
    os.makedirs(megahit_subset, exist_ok=True)
    os.makedirs(spades_subset, exist_ok=True)
    if not args.no_subset:
        run_command(f"zcat {args.megahit_contigs} | seqtk sample -s100 - {args.num_samples} > {megahit_subset}/megahit_subset.fasta")
        run_command(f"zcat {args.spades_contigs} | seqtk sample -s100 - {args.num_samples} > {spades_subset}/spades_subset.fasta")

    # Run MaxBin2
    megahit_bins = os.path.join(args.output_dir, 'megahit_bins')
    spades_bins = os.path.join(args.output_dir, 'spades_bins')
    os.makedirs(megahit_bins, exist_ok=True)
    os.makedirs(spades_bins, exist_ok=True)
    if not args.no_maxbin:
        run_command(f"metabat2 -i {megahit_subset}/megahit_subset.fasta -o {megahit_bins}/megahit_bins")
        run_command(f"metabat2 -i {spades_subset}/spades_subset.fasta -o {spades_bins}/spades_bins")

    # Run CheckM
    if not args.no_checkm:
        megahit_checkm = os.path.join(args.output_dir, 'megahit_checkm')
        spades_checkm = os.path.join(args.output_dir, 'spades_checkm')
        run_command(f"checkm lineage_wf {megahit_bins} {megahit_checkm} -x fa -t 10")
        run_command(f"checkm lineage_wf {spades_bins} {spades_checkm} -x fa -t 10")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Subset contigs and run MaxBin2 and CheckM")
    parser.add_argument("megahit_contigs", help="Path to Megahit contigs (gzipped)")
    parser.add_argument("spades_contigs", help="Path to Spades contigs (gzipped)")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for bins and CheckM results")
    parser.add_argument("--num_samples", type=int, default=1000, help="Number of samples to subset (default=1000)")
    parser.add_argument("--no-subset", action="store_true", help="Skip the subset step")
    parser.add_argument("--no-maxbin", action="store_true", help="Skip the MaxBin2 step")
    parser.add_argument("--no-checkm", action="store_true", help="Skip the CheckM step")

    args = parser.parse_args()
    main(args)
