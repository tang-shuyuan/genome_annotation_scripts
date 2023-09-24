import pysam
import datetime
import sys
import os
import multiprocessing
import pyBigWig
import argparse
parser=argparse.ArgumentParser(description="this is a python script used to merge bam files become a bigwig file")
parser.add_argument('input_folder_or_files',nargs='+', help='input bam files or a folder containing BAM files ')
parser.add_argument('-o', '--output_file', default='output.bigwig', help='Output a bigwig file')
parser.add_argument('-t','--cpu' , type=int, default=1, help='default=1 The number of CPUs used during code execution')
args = parser.parse_args()

input_bam_files=[]
for item in args.input_folder_or_files:
    if os.path.isdir(item):
        input_bam_files += [os.path.join(item, filename) for filename in os.listdir(item) if filename.endswith('.bam')]
    elif os.path.isfile(item) :
        input_bam_files.append(item)
    else:
        print(f"'{item}' is not a valid file or folder")
if not input_bam_files:
    print("No BAM files found")
    sys.exit()
output_file = args.output_file

current_time = datetime.datetime.now()
print("current_time：", current_time)
print("the following bamfile will be merge to a bigwig file:")
for bam in input_bam_files:
    print(bam)
print(f"used cpu number is(are) {args.cpu}")

def get_chr_imformations(bam_files):
    chromosomes = {}
    with pysam.AlignmentFile(bam_files[0], "rb") as bam:
        for contig in bam.header['SQ']:
            chromosomes[contig['SN']] = contig['LN']
    return chromosomes

def calculate_coverage(input_bam, chr_name, chr_length):
    coverage = [0] * chr_length
    with pysam.AlignmentFile(input_bam, "rb") as bam:
        for read in bam.fetch(chr_name):
            positions = read.get_reference_positions(full_length=True)
            positions = [pos for pos in positions if pos is not None]
            for position in positions:
                try:
                    coverage[position] += 1
                except IndexError:
                    print(f"{position} out of index, don't worry, ignore this")
    return coverage
def merged_bamfiles_coverage(bamfiles, chr_name, chr_length):
    merged_coverage = [0] * chr_length
    pool = multiprocessing.Pool(processes=args.cpu)
    results = []
    for bamfile in bamfiles:
        result = pool.apply_async(calculate_coverage, (bamfile, chr_name, chr_length))
        results.append(result)
    for result in results:
        coverage = result.get()
        merged_coverage = [sum(x) for x in zip(merged_coverage, coverage)]
    pool.close()
    pool.join()
    return merged_coverage

def process_chr_and_write_bigwig(bamfiles,output_file,chr_imformations):
    bw = pyBigWig.open(output_file, 'w')
    bw.addHeader(list(chr_imformations.items()))
    for chr_name,chr_length in chr_imformations.items():
        coverage = merged_bamfiles_coverage(bamfiles,chr_name,chr_length)
        i = 0
        while i < chr_length:
            a = coverage[i]
            b = i + 1
            while b < chr_length and coverage[b] == a:
                b += 1
            else:
                bw.addEntries([chr_name], [i], ends=[b], values=[round(float(a), 1)])
                i = b
    bw.close()
chr_imformations=get_chr_imformations(input_bam_files)
process_chr_and_write_bigwig(input_bam_files, output_file,chr_imformations)
current_time = datetime.datetime.now()
print("current_time:", current_time)