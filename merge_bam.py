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
parser.add_argument('-r','--interval' , type=int, default=5000, help='default=5000 provide a number as the distance to split non gene region ')
parser.add_argument('-o', '--output_dir', default='output', help='Output directory')
parser.add_argument('-s','--sort_and_index',action='store_true',help='Whether to sort and index Bam files')
parser.add_argument('-t','--cpu' , type=int, default=4, help='default=4 The number of CPUs used during code execution')
parser.add_argument('-m','--max_intron_length',type=int,default=100000,help='default=100kb,Filter reads containing abnormal intron length')
args = parser.parse_args()
bam_file_paths=[]
for item in args.input_folder_or_files:
    if os.path.isdir(item):
        bam_file_paths += [os.path.join(item, filename) for filename in os.listdir(item) if filename.endswith('.bam')]
    elif os.path.isfile(item) :
        bam_file_paths.append(item)
    else:
        print(f"'{item}' is not a valid file or folder")

if len(bam_file_paths) <= 1:
    print("at least two bam files need to be merged")
    sys.exit()
print("The following bam files will be merged:")
for bam in bam_file_paths:
    print(bam)

if args.interval <0:
    print("The -r parameter needs to be set correctly")
    sys.exit()

input_gff_file=args.gff_file
interval=args.interval
process_num=args.cpu
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)
max_intron_length=args.max_intron_length   #max intron length


start_time=time.time()
current_time = datetime.datetime.now()
print("current_time：", current_time)

print(f"the used cpu number is(are) {process_num}")
def sort_and_index_bamfile(bamfile, output_dir):
    filename = os.path.basename(bamfile)
    sorted_bam = os.path.join(output_dir, f"sorted_{filename}")
    pysam.sort("-o", sorted_bam, bamfile)
    print(f"{filename} sorted successfully")
    pysam.index(sorted_bam)
    print(f"{filename} indexed successfully")

def parallel_sort_and_index_bamfiles(bamfiles, output_dir, num_processes):
    pool = multiprocessing.Pool(processes=num_processes)
    # Use pool.starmap to pass multiple arguments to sort_and_index_bamfile
    pool.starmap(sort_and_index_bamfile, [(bamfile, output_dir) for bamfile in bamfiles])
    pool.close()
    pool.join()

if args.sort_and_index:
    print("Sorting and indexing BAM files in parallel")
    # num_processes = 4  # You can adjust the number of processes based on your system's capabilities
    parallel_sort_and_index_bamfiles(bam_file_paths, output_dir, process_num)
    bam_file_paths = glob.glob(output_dir + "/*bam")
    print("sorting and indexing complete")
    current_time = datetime.datetime.now()
    print("current_time：", current_time)


# #####Get chromosome name and length information
def get_chrom_lengths(bam_file):
    chrom_lengths = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for contig in bam.header["SQ"]:
            chrom_lengths[contig["SN"]] = contig["LN"]
    bam.close()
    return chrom_lengths
# #run get_chrom_lengths，get chr，length
chr_len=get_chrom_lengths(bam_file_paths[0])


# # Read the GFF file and get the gene region
def split_gff_to_bed(gff_file,chrom_lengths):
    #get gene region
    gene_regions = []
    with open(gff_file, "r") as file:
        for line in file:
            if not line.startswith("#") and line.strip():
                columns = line.strip().split("\t")
                feature_type = columns[2].lower()
                if feature_type in ["gene", "mrna", "transcript"]:
                    chromosome = columns[0]
                    start = int(columns[3])
                    end = int(columns[4])
                    gene_regions.append((chromosome, start, end))
    # mergeing overlapping gene regions
    merged_regions = []
    gene_regions.sort(key=lambda x: (x[0], x[1]))  # Sorted by chromosome and gene starting position
    for region in gene_regions:
        if merged_regions:
            last_region = merged_regions[-1]
            if region[0] == last_region[0] and region[1] <= last_region[2]:
                last_region = (last_region[0], last_region[1], max(region[2], last_region[2]))
                merged_regions[-1] = last_region
                continue
        merged_regions.append(region)


# # ####get nongene_region
    non_gene_regions = []
    merged_regions.sort(key=lambda x: (x[0], x[1]))
    for chrom,lenght in chrom_lengths.items():
        start=0
        for chr, gene_start, gene_end in merged_regions:
            if chr == chrom and start > gene_start:
                print("error number")
            while chr==chrom and start<gene_start:
                if gene_start-start>interval:
                        non_gene_regions.append((chrom, start, start+interval))
                        start = start+interval
                elif gene_start-start<=interval:
                    non_gene_regions.append((chrom, start, gene_start))
                    start=gene_end
        non_gene_regions.append((chrom,start,lenght))
    all_regions = merged_regions + non_gene_regions
    with open(os.path.join(output_dir, "region.txt"),"w") as f:
        merged_regions = [(region + ("gene_region",)) for region in merged_regions]
        non_gene_regions = [(region + ("non_gene_region",)) for region in non_gene_regions]
        list_region = sorted(merged_regions + non_gene_regions, key=lambda x: (x[0], x[1]))
        f.writelines([f"{content}\n" for content in list_region])
    return merged_regions,all_regions

gene_region_list,unsort_all_region_list=split_gff_to_bed(input_gff_file,chr_len)
all_region_list=sorted(unsort_all_region_list,key=lambda x: (x[0], x[1]))
###Check whether the chromosome names correspond
store_set=set()
for tuple in gene_region_list:
    store_set.add(tuple[0])
for chr_name in chr_len.keys():
    if chr_name not in store_set:
        print(f"warning, {chr_name} cannot be found in the gff file,this chromosome will not be split region")
print("split region succssfully")


# ###calculate the  depth of each region
def calculate_coverage(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    coverage_dict = {}
    for region in all_region_list:
        chrom, start, end = region
        pileup = bam.fetch(chrom, start, end)
        depth = sum([1 for _ in pileup])
        coverage_dict[region] = depth
    bam.close()
    return coverage_dict

def process_bam_file(bam):
    coverage = calculate_coverage(bam)
    print("collect "+ bam + " depth ifomation  finished")
    return bam, coverage
def store_bam_region_coverage_parallel(bam_files):
    now_time = time.time()
    print("Start calculating the region depth of each bam file in parallel")
    print(f"by far ,spent time {now_time - start_time} seconds")
    coverage_dicts = {}
    with multiprocessing.Pool() as pool:
        results = [pool.apply_async(process_bam_file, (bam,)) for bam in bam_files]
        for result in results:
            bam, coverage = result.get()
            coverage_dicts[bam] = coverage
    end_time = time.time()
    total_time = end_time - start_time
    print(f"The depth information is collected and the total time for parallel computation is {total_time} seconds")
    return coverage_dicts
big_dict=store_bam_region_coverage_parallel(bam_file_paths)

def merge_regions_to_bam(bam_files,coverage_dicts,output_region):
    create_header_bam=pysam.AlignmentFile(bam_files[0],"rb")
    merged_bam = pysam.AlignmentFile(os.path.join(output_dir, "merged.bam"), "wb", header=create_header_bam.header)
    def comprise_max_bam(region):
        max_coverage = 0
        max_bam = None
        for bam in coverage_dicts:
            coverage=coverage_dicts[bam][region]
            if coverage>=max_coverage:
                max_coverage=coverage
                max_bam = bam
        return max_bam
    region_number=len(all_region_list)
    store_dict = {}
    i = 0
    while i < region_number:
        initial_max_bam=comprise_max_bam(all_region_list[i])
        b=i+1
        continuous_region=all_region_list[i]
        while comprise_max_bam(all_region_list[b])==initial_max_bam and all_region_list[b][0]==all_region_list[i][0]:
            continuous_region=(all_region_list[i][0],all_region_list[i][1],all_region_list[b][2])
            b+=1
            if b==region_number:
                break
        store_dict[continuous_region]=initial_max_bam
        i=b
        if i == region_number-1:
            store_dict[all_region_list[-1]]=comprise_max_bam(all_region_list[-1])
            i=i+1
    with open(output_region,"w") as f:
        for key,values in store_dict.items():
            print((key,values),file=f)
    for bam in bam_files:
        region_merged_bam = pysam.AlignmentFile(bam, "rb")
        for region,bamname in store_dict.items():
            if bam == bamname:
                for read in region_merged_bam.fetch(*region):
                    ###limit intron length
                    contain_big_intron = False
                    if read.cigartuples is not None:
                        contain_big_intron = any(op == 3 and length >= max_intron_length for op, length in read.cigartuples)
                    if not contain_big_intron:
                        merged_bam.write(read)
        region_merged_bam.close()
    merged_bam.close()
    print("completed mergeing bam files")
merge_regions_to_bam(bam_file_paths,big_dict,os.path.join(output_dir,"merged_region.txt"))

###进行最后排序和建索引
print("Start sorting and indexing the merged bam files")
def sort_and_index_merged_bam(input_bam,sorted_bam,index_file):
    # Sort BAM file
    pysam.sort("-o", sorted_bam, input_bam)
    # Create index for sorted BAM file
    pysam.index(sorted_bam, index_file)
sort_and_index_merged_bam(os.path.join(output_dir, "merged.bam"), os.path.join(output_dir,"merged_sorted.bam"),os.path.join(output_dir,"merged_sorted.bam.bai"))
print("Complete sorting and indexing ")
###建立索引和排序好后删除merged.bam
# try:
#     os.remove(os.path.join(output_dir, "merged.bam"))
#     print("merged.bam deleted")
# except FileNotFoundError:
#     print("file not exist")
# except Exception as e:
#     print(f"An error occurred while deleting the file：{e}")

end_time = time.time()
duration = end_time - start_time
print(f"scipt runtime {duration} seconds")
current_time = datetime.datetime.now()
print("current_time:", current_time)
