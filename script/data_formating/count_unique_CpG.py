import os
import argparse
import gzip

def count_unique_cpg_sites(bed_folder, output_file):
    cpg_sites = set()
    bed_files = [f for f in os.listdir(bed_folder) if f.endswith(".bed") or f.endswith(".bed.gz")]

    for bed_file in bed_files:
        file_path = os.path.join(bed_folder, bed_file)
        
        if bed_file.endswith(".gz"):
            with gzip.open(file_path, 'rt') as file:
                for line in file:
                    if line.startswith('#') or line.strip() == "":
                        continue
                    cols = line.strip().split('\t')
                    chrom = cols[0]
                    start = cols[1]
                    end = cols[2]
                    cpg_sites.add((chrom, start, end))
        else:
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith('#') or line.strip() == "":
                        continue
                    cols = line.strip().split('\t')
                    chrom = cols[0]
                    start = cols[1]
                    end = cols[2]
                    cpg_sites.add((chrom, start, end))
    return len(cpg_sites)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compte unique CpG sites and create a list")
    parser.add_argument("bed", help="Dossier contenant les fichiers BED")
    args = parser.parse_args()

    unique_cpg_count = count_unique_cpg_sites(args.bed)
    print(f"Nombre de sites CpG uniques : {unique_cpg_count}")