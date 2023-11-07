#!/usr/bin/env python3

# cannot get the conda or homebrew quast installs to work, which seems weird. installed quast from source.
# conda install -y -c bioconda pilon
# conda install -y -c bioconda hypo
# conda install -y -c bioconda bbmap ### for tadpole
import argparse
from tqdm import tqdm
from Bio import SeqIO
import subprocess
import random
import os

def write_subset(input_file, output_file, subset_size):
    records = list(SeqIO.parse(input_file, "fasta"))
    random.shuffle(records)
    subset = random.sample(records, min(subset_size, len(records)))
    SeqIO.write(subset, output_file, "fasta")

def run_tool(tool, input_file, output_file):
    if tool == "pilon":
        cmd = f"pilon --genome {input_file} --output {output_file}"
    elif tool == "hypo":
        cmd = f"hypo -i {input_file} -o {output_file}"
    elif tool == "tadpole":
        cmd = f"tadpole.sh in={input_file} out={output_file}"
    else:
        print(f"Unknown tool: {tool}")
        return
    subprocess.run(cmd, shell=True)

def run_quast(corrected_file, original_file, quast_output_dir):
    cmd = f"quast -o {quast_output_dir} {corrected_file} -r {original_file}"
    subprocess.run(cmd, shell=True)

def parse_quast_report(quast_dir):
    report_path = os.path.join(quast_dir, 'report.tsv')
    metrics = {}
    try:
        with open(report_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:  # Skip the header
                parts = line.strip().split('\t')
                metric, value = parts[0], parts[1]
                metrics[metric] = value
    except FileNotFoundError:
        print(f"QUAST report not found in {quast_dir}")
    return metrics

def main():
    parser = argparse.ArgumentParser(description="Test assembly correction tools.")
    parser.add_argument("-m", "--megahit_contigs", required=True, help="MEGAHIT contigs file")
    parser.add_argument("-s", "--spades_contigs", required=True, help="SPAdes contigs file")
    parser.add_argument("-n", "--subset_size", type=int, default=10000, help="Number of contigs for subset")
    args = parser.parse_args()

    tools = ["pilon", "hypo", "tadpole"]
    quast_results = {}

    for assembler, contigs in tqdm({"megahit": args.megahit_contigs, "spades": args.spades_contigs}.items()):
        subset_file = f"{assembler}_subset.fasta"
        write_subset(contigs, subset_file, args.subset_size)

        for tool in tools:
            output_file = f"{assembler}_{tool}_corrected.fasta"
            run_tool(tool, subset_file, output_file)

            # Run QUAST and parse results
            quast_dir = f"{assembler}_{tool}_quast_output"
            run_quast(output_file, subset_file, quast_dir)
            quast_metrics = parse_quast_report(quast_dir)
            quast_results[f"{assembler}_{tool}"] = quast_metrics

    # Output QUAST comparison
    print("QUAST Metrics Comparison:")
    for key, metrics in quast_results.items():
        print(f"{key}: {metrics}")
if __name__ == "__main__":
    main()
