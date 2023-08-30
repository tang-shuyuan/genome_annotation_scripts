import sys
import os
import glob
import multiprocessing
import argparse
import pysam
import time
import datetime
parser = argparse.ArgumentParser(description="this is a python script used to merge bam files")
parser.add_argument('input_folder_or_files',nargs='+', help='input a folder containing BAM files or bam files')
parser.add_argument('gff_file', help='Input GFF file')
parser.add_argument('-c','--coverage' , type=float, default=70, help='default=70 The coverage rate threshold for a gene to be considered as expressed(Number size must be between 0 and 100')
parser.add_argument('-d','--depth' , type=float, default=10, help='default=10 The coverage depth threshold for a gene to be considered as expressed(non-negative)')
parser.add_argument('-o', '--output_file', default='output_file', help='Output directory')
parser.add_argument('-t','--cpu' , type=int, default=1, help='default=1 The number of CPUs used during code execution')
args = parser.parse_args()
bam_file_paths=[]
for item in args.input_folder_or_files:
    if os.path.isdir(item):
        bam_file_paths += [os.path.join(item, filename) for filename in os.listdir(item) if filename.endswith('.bam')]
    elif os.path.isfile(item) :
        bam_file_paths.append(item)
    else:
        print(f"'{item}' is not a valid file or folder")

if not bam_file_paths:
    print("No BAM files found")
    sys.exit()

if args.coverage <0 or args.coverage>100:
    print("the coverage must be between 0 and 100")
    sys.exit()

if args.depth <0 :
    print("The -d parameter needs to be set correctly")
    sys.exit()

input_gff_file=args.gff_file
process_num=args.cpu
output_file = args.output_file


start_time=time.time()
current_time = datetime.datetime.now()
print("current_time：", current_time)
print(f"the used cpu number is(are) {process_num}")

def calculate_all_genome_coverage_depth(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    total_bases = 0
    covered_bases = 0
    covered_depth=0
    for ref in bam.references:
        total_bases += bam.lengths[bam.references.index(ref)]
    for pileupcolumn in bam.pileup():
        if pileupcolumn.nsegments > 0:
            covered_depth+=pileupcolumn.nsegments
            covered_bases += 1
    bam.close()
    coverage = round(covered_bases / total_bases*100,2)
    depth= round(covered_depth / total_bases,2)
    return coverage,depth

def calculate_gene_coverage_depth(bamfile,gff_file):
    bam = pysam.AlignmentFile(bamfile, "rb")
    expressed_gene_number=0
    not_expressed_gene_number=0
    with open(gff_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                columns = line.strip().split("\t")
                feature_type = columns[2]
                if feature_type == "gene":
                    chromosome = columns[0]
                    start = int(columns[3])
                    end = int(columns[4])
                    gene_length= end-start+1
                    covered_bases = 0
                    covered_depth = 0
                    for pileupcolumn in bam.pileup(chromosome, start, end):
                        if pileupcolumn.reference_pos >= start and pileupcolumn.reference_pos <= end:
                            if pileupcolumn.n > 0:
                                covered_depth+=pileupcolumn.n
                                covered_bases += 1
                    coverage=round(covered_bases/gene_length*100,2)
                    depth=round(covered_depth/gene_length,2)
                    if depth>args.depth and coverage>args.coverage:
                        expressed_gene_number+=1
                    else:
                        not_expressed_gene_number+=1
        total_gene_number=expressed_gene_number+not_expressed_gene_number
    bam.close()
    return total_gene_number,expressed_gene_number
def process_bam_file(bam_file, input_gff_file):
    coverage, depth = calculate_all_genome_coverage_depth(bam_file)
    total_gene_number, expressed_gene_number = calculate_gene_coverage_depth(bam_file, input_gff_file)
    return bam_file, coverage, depth, total_gene_number, expressed_gene_number
current_time = datetime.datetime.now()
print("current_time:", current_time)
print("start calculate whole genome coverage and gene expression levels for each bamfile")
if __name__ == '__main__':
    with open(output_file, 'w') as f:
        sys.stdout = f
        pool = multiprocessing.Pool(processes=process_num)
        results = []
        for bam_file in bam_file_paths:
            results.append(pool.apply_async(process_bam_file, (bam_file, input_gff_file)))
        pool.close()
        pool.join()
        for result in results:
            bam_file, coverage, depth, total_gene_number, expressed_gene_number = result.get()
            print(f"The coverage(ratio) of {bam_file} is {coverage}%")
            print(f"The depth of {bam_file}' is {depth}")
            print(f"{bam_file}'s total gene number are {total_gene_number},expressed gene number are {expressed_gene_number}")
        sys.stdout = sys.__stdout__

end_time = time.time()
duration = end_time - start_time
print(f"scipt runtime {duration} seconds")
current_time = datetime.datetime.now()
print("current_time:", current_time)
