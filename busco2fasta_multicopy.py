#!/usr/bin/env python3
import sys
import argparse
import os
import shutil

def get_list_of_busco_dirs(results_dir):
    # get list of busco_dirs
    busco_dir_list = []
    for busco_dir in os.listdir(results_dir):
        busco_dir_list.append(busco_dir)
    # print output
    print("[STATUS] Found " + str(len(busco_dir_list)) + " BUSCO runs in input directory ('" + results_dir.replace("/", "") + "'):")
    for busco_dir in sorted(busco_dir_list):
        print("[STATUS] \t- " + busco_dir)
    # return list
    return busco_dir_list

def count_single_copy_buscos(results_dir, busco_dir_list, suffix):
    scb_count_dict = {}
    for busco_dir in busco_dir_list:
        for item in os.listdir(results_dir + busco_dir):
            if item.startswith("run_"):
                lineage = item.replace("run_", "")
                # Define the directories to check
                directories = {
                    "single_copy": results_dir + busco_dir + "/run_" + lineage + "/busco_sequences/single_copy_busco_sequences/",
                    "multi_copy": results_dir + busco_dir + "/run_" + lineage + "/busco_sequences/multi_copy_busco_sequences/"
                }
                # Iterate over the directories
                for copy_type, directory in directories.items():
                    try:
                        scb_file_list = os.listdir(directory)
                        for scb_file in scb_file_list:
                            if scb_file.endswith(suffix):
                                scbID = scb_file.replace(suffix, "")
                                try:
                                    scb_count_dict[scbID] += 1
                                except KeyError:
                                    scb_count_dict[scbID] = 1
                                # Print the basename of the file being processed
                                print(f"Processing file: {os.path.basename(scb_file)}")
                                # Process the file and append the appropriate suffix
                                process_file(directory, scb_file, copy_type)
                    except FileNotFoundError:
                        continue
    return scb_count_dict, lineage

def process_file(directory, filename, copy_type):
    suffix = "sc" if copy_type == "single_copy" else "mc"
    filepath = os.path.join(directory, filename)
    
    # Read the content of the file
    with open(filepath, 'r') as infile:
        lines = infile.readlines()
    
    # Modify the headers
    modified_lines = []
    for line in lines:
        if line.startswith(">"):
            modified_lines.append(line.strip() + f"{suffix}\n")
        else:
            modified_lines.append(line)
    
    return modified_lines

def get_usable_buscos(scb_count_dict, proportion, busco_dir_list):
    n_taxa = len(busco_dir_list)
    min_taxon_count = int(round(n_taxa * proportion, 0))
    usable_scb_list = []
    print("[STATUS] Selecting BUSCOs present (single-copy OR multi_copy) in ≥ " + str(min_taxon_count) + " taxa...")
    for scbID, count in scb_count_dict.items():
        if count >= min_taxon_count:
            usable_scb_list.append(scbID)
    return usable_scb_list

def parse_fasta(fasta_file, prefix, suffix):
    with open(fasta_file) as fasta:
        fasta_dict = {}
        for line in fasta:
            if line.startswith(">"):
                header = ">" + prefix + "." + line.rstrip("\n").replace(">", "").split(" ")[0] + suffix
                fasta_dict[header] = ''
            else:
                fasta_dict[header] += line.rstrip("\n")
    return fasta_dict

def create_output_fastas(results_dir, busco_dir_list, usable_scb_list, suffix, outdir, lineage):
    print("[STATUS] Writing " + str(len(usable_scb_list)) + " " + seqtype + " FASTA files to output directory ('" + (outdir) + "'). This may take a moment...") 
    total, count = len(usable_scb_list), 0
    for scbID in usable_scb_list:
        output_fasta_list = []
        for busco_dir in busco_dir_list:
            for copy_type in ["single_copy", "multi_copy"]:
                path = results_dir + busco_dir + "/run_" + lineage + f"/busco_sequences/{copy_type}_busco_sequences/"
                try:
                    fasta_dict = parse_fasta(path + scbID + suffix, busco_dir, "_sc" if copy_type == "single_copy" else "_mc")
                    output_fasta_list.append(fasta_dict)
                except FileNotFoundError:
                    pass
        with open(outdir + "/" + scbID + suffix, 'a') as outfile:  # Open in append mode
            for fasta_dict in output_fasta_list:
                for header, sequence in fasta_dict.items():
                    outfile.write(header + "\n")
                    outfile.write(sequence + "\n")
        print("[STATUS] \t" + str(count) + "/" + str(total) + " files written.", end='\r')
        count += 1
    print("[STATUS] \t" + str(count) + "/" + str(total) + " files written.")

if __name__ == "__main__":
    SCRIPT = "busco2fasta.py"
    ### argument set up 
    # set up
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--busco_dir", type=str, help = "directory containing a set BUSCO results directories (e.g. one per taxon)", required=True)
    parser.add_argument("-o", "--outdir", type=str, help = "output directory for FASTAs (default: b2f_output)", default="b2f_output")
    parser.add_argument("-s", "--seqtype", type=str, help = "your chosen sequence type (default: protein)", choices={"protein", "nucleotide"}, default="protein")
    parser.add_argument("-p", "--proportion", type=float, help = "proportion of taxa required for a given BUSCO to be output as FASTA (default: 1.0)", default=1.0)
    args = parser.parse_args()
    results_dir = args.busco_dir
    seqtype = args.seqtype
    proportion = args.proportion
    outdir = args.outdir
    # add slash to results dir if necessary
    if not results_dir.endswith("/"):
        results_dir = results_dir + "/"
    # get suffix based on seqtype    
    if seqtype == 'protein':
        suffix = ".faa"
    else:
        suffix = ".fna"
    # proportion
    print("busco2fasta.py parameters:")
    print("\tinput_dir: " + results_dir.replace("/", ""))
    print("\toutput_dir: " + outdir)
    print("\tseqtype: " + seqtype)
    print("\tproportion: " + str(proportion))
    print("")
    if os.path.exists(outdir):
        print("[WARNING] Removing existing output directory ('" + outdir + "') and its contents...")
        shutil.rmtree(outdir)
    os.mkdir(outdir)
    ### run functions    
    busco_dir_list = get_list_of_busco_dirs(results_dir)
    scb_count_dict, lineage = count_single_copy_buscos(results_dir, busco_dir_list, suffix)
    usable_scb_list = get_usable_buscos(scb_count_dict, proportion, busco_dir_list)
    create_output_fastas(results_dir, busco_dir_list, usable_scb_list, suffix, outdir, lineage)
    print("[STATUS] busco2fasta.py completed successfully.")
