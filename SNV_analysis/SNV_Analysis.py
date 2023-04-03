import sys
import subprocess
import pickle
import os
from collections import defaultdict
import time

def occurence_filter(current_dict, threshold):
    """
    current_dict: previous threshold dictionary
    threshold: proportion to filter at (0-1.0)
    
    return a filtered version of the dictionary 
    """
    new_dict = dict()
    for key in current_dict:
        if int(current_dict[key][key[2]]) > int(threshold * (current_dict[key]['A'] + 
                                                             current_dict[key]['T'] + 
                                                             current_dict[key]['G'] + 
                                                             current_dict[key]['C'] + 
                                                             current_dict[key]['Deletions'])):
            new_dict[key] = current_dict[key]
    return new_dict


def make_gvf_set(gvf_file_directory, min_snv_dict):
    """
    gvf_file_directory: directory with Ensembl gvf files
    min_snv_dict: For the current sample, the minimum threshold (30%) dictionary
    
    return a set of all known variants
    """
    gvf_set = set()
    snv = set([(key[0], key[1]) for key in min_snv_dict])
    for file in os.listdir(gvf_file_directory):
        chrom = file.split("-")[-1].split(".")[0]
        if "homo_sapien" not in file:
            continue
        else:
            for line in open(os.path.join(gvf_file_directory, file)):
                if "SNV" not in line:
                    continue
                split_line = line.split()
                pos = int(split_line[3])
                if (chrom, pos) not in snv:
                    continue
                info = split_line[-1].split(";")
                var_seq = ''
                for element in info:
                    if 'Variant_seq' in element:
                        var_seq = element.split("=")[-1]
                        if len(var_seq.split(',')) > 1:
                            for t in var_seq.split(','):
                                gvf_set.add((chrom, pos, t))
                        else:
                            gvf_set.add((chrom, pos, var_seq))
    return gvf_set


def snvs_in_gvf(gvf, snv_dict):
    """
    gvf: set of known Ensembl variants
    snv_dict: current snv dict to bin
    
    returns a tuple with the count of gvf filtered, the dictionary with gvf sites filtered out
    and then the dictionary with the gvf sites
    """
    accounted_for = 0
    filtered_dict = {}
    gvf_dict = {}
    for key in list(snv_dict.keys()):
        if key in gvf:
            accounted_for += 1
            gvf_dict[key] = snv_dict[key]
        else:
            filtered_dict[key] = snv_dict[key]
    return (accounted_for, filtered_dict, gvf_dict)

def get_ref_seq(path_to_ref, ref_section, start_on_ref, end_on_ref):
    """
    path_to_ref: path to reference sequence
    ref_section: contig (chromosome for complete annotated references)
    start_on_ref: start position to extract
    end_on_ref: ending position to extract
    
    return string of reference sequence
    """
    ref_seq = subprocess.run(['samtools', 
                              'faidx', 
                              f"{path_to_ref}", 
                              f"{ref_section}:{start_on_ref}-{end_on_ref}"],
                             stdout=subprocess.PIPE).stdout.decode('utf-8')
    ref_seq = ref_seq.split('\n')[1:]
    ref_seq = "".join(ref_seq)
    return ref_seq


def count_prob_kmers(snv_dict, prob_kmer_set, p_to_ref):
    """
    snv_dict: current dictionary of snvs
    prob_kmer_set: low confidence kmers
    p_to_ref: path to the reference file
    
    returns a tuple with the count of low confidence kmers filtered, 
    the dictionary with low confidence kmer sites filtered out
    and then the dictionary with the low confidence kmer sites
    """
    count = 0
    prob_kmers = {}
    filtered_dict = {}
    for key in snv_dict:
        chrom = key[0]
        pos = int(key[1])
        kmer = get_ref_seq(p_to_ref, chrom, pos - 4, pos + 4)
        if kmer in prob_kmer_set:
            count += 1
            prob_kmers[key] = snv_dict[key]
        else:
            filtered_dict[key] = snv_dict[key]
    return (count, filtered_dict, prob_kmers)

def snv_dict_value_to_tab_string(snv_dict_value):
    """
    snv_dict_value:Values for SNV_dict
    
    return tab delimited form of SNV_dict value
    """
    s = ""
    for key in snv_dict_value:
        if key == 'Reads':
            continue
        s = s + f"{snv_dict_value[key]}\t"
    s = s.rstrip()
    return s  


def variant_bin_and_write_to_file(snv_dict, file_id, gvf_set, kmer_set, path_to_ref, out_directory_path):
    """
    snv_dict: the snv dict to bin variants from
    file_id: tag for the output (usually the occurence threshold ex. v30)
    gvf_set: set of Ensembl known variant locations
    kmer_set: set of Low confidence kmers
    path_to_ref: path to the reference genome
    out_directory_path: Path to store output files
    
    return N/A
    """
    #Filter based on known variants
    gvf_filter_results = snvs_in_gvf(gvf_set, snv_dict)
    gvf_dict = gvf_filter_results[2]
    gvf_filtered_snv_dict = gvf_filter_results[1]
    gvf_count = gvf_filter_results[0]
    
    #Filter based on problematic Kmers
    prob_kmers_results = count_prob_kmers(gvf_filtered_snv_dict, kmer_set, path_to_ref)
    prob_kmer_dict = prob_kmers_results[2]
    final_snv_dict = prob_kmers_results[1]
    prob_kmer_count = prob_kmers_results[0]
    
    #Write Known variations to file
    fh = open(f"{out_directory_path}/{file_id}_known_variants.tsv", 'w')
    fh.write("Chrom\tPosition\tVariant\tReference\tA_count\tC_count\tT_count\tG_count\tDel_count\tInsert_count\n")
    for key in gvf_dict:
        line = f"{key[0]}\t{key[1]}\t{key[2]}\t{snv_dict_value_to_tab_string(gvf_dict[key])}\n"
        fh.write(line)
    fh.close()
    
    #Write Problem kmer variations to file
    fh = open(f"{out_directory_path}/{file_id}_problem_kmer_variants.tsv", 'w')
    fh.write("Chrom\tPosition\tVariant\tReference\tA_count\tC_count\tT_count\tG_count\tDel_count\tInsert_count\n")
    for key in prob_kmer_dict:
        line = f"{key[0]}\t{key[1]}\t{key[2]}\t{snv_dict_value_to_tab_string(prob_kmer_dict[key])}\n"
        fh.write(line)
    fh.close()
    
    #Write novel variations to file
    fh = open(f"{out_directory_path}/{file_id}_novel_IVT_variants.tsv", 'w')
    fh.write("Chrom\tPosition\tVariant\tReference\tA_count\tC_count\tT_count\tG_count\tDel_count\tInsert_count\n")
    for key in final_snv_dict:
        line = f"{key[0]}\t{key[1]}\t{key[2]}\t{snv_dict_value_to_tab_string(final_snv_dict[key])}\n"
        fh.write(line)
    fh.close()

    
def main(pysam_file_path, ref_file_path, gvf_file_directory, out_directory_path):
    """
    pysam_file_path: pysamstats file w/ columns: chrom,pos,ref,reads_all,
                                                 deletions,insertions,A,C,
                                                 T,G,N                                  
    ref_file_path: path to the reference file the alignments were originally
                   performed on
    gvf_file_directory: Directory containing the Ensembl gvf files
    out_directory_path: Output directory for variants
    
    return N/A                                            
    """
    
    
    options = ["A", "C", "T", "G"]
    snv_dict30 = dict()
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 
              'chr5', 'chr6', 'chr7', 'chr8', 
              'chr9','chr10', 'chr11', 'chr12',
              'chr13', 'chr14', 'chr15', 'chr16', 
              'chr17', 'chr18', 'chr19', 'chr20',
              'chr21','chr22', 'chrX', 'chrM']
    
    print("Reading pysamstats in")
    for line in open(pysam_file_path, 'r'):
        
        #Skip the header line
        if 'reads_all' in line:
            continue
        split_line = line.split()
        chrom = split_line[0]

        #Omit Mitochondrial and Y chromosome
        if chrom not in chroms:
            continue
        
        pos = int(split_line[1])
#        print(f"{chrom} {pos}")
        ref = split_line[2]
        Deletions = int(split_line[4])
        Insertions = int(split_line[5])
        
        option_dict = {"A": int(split_line[6]), 
                       "C": int(split_line[7]),
                       "T": int(split_line[8]),
                       "G": int(split_line[9])
                      }
        
        # Exclude Insertions in establshing number of reads at position
        reads = (option_dict["A"] + 
                 option_dict["T"] + 
                 option_dict["C"] + 
                 option_dict["G"] + 
                 Deletions)
        
        # Only consider positions with 10 reads
        if reads < 10:
            continue
        
        for option in ["A", "T", "C", "G"]:
            
            #Only consider options that are possible SNVs (not the ref)
            if option == ref:
                continue
            if option_dict[option] >= reads * 0.3:
                snv_dict30[(chrom, pos, option)] = {"Ref": ref,
                                                    "A": option_dict["A"],
                                                    "C": option_dict["C"],
                                                    "T": option_dict["T"],
                                                    "G": option_dict["G"],
                                                    "Deletions": Deletions,
                                                    "Insertions": Insertions}
        
    #Make the gvf set
    print("gvf")
    gvf_set = make_gvf_set(gvf_file_directory, snv_dict30)
        
    #Load the precomputed low conf kmer set
    fh = open('low_conf_kmer', 'rb')
    low_conf_kmer_set = pickle.load(fh)
    fh.close()
    
    #Compute the occurence thresholds
    print("Dicts")
    snv_dict40 = occurence_filter(snv_dict30, 0.4)
    snv_dict50 = occurence_filter(snv_dict40, 0.5)
    snv_dict60 = occurence_filter(snv_dict50, 0.6)
    snv_dict70 = occurence_filter(snv_dict60, 0.7)
    snv_dict80 = occurence_filter(snv_dict70, 0.8)
    snv_dict95 = occurence_filter(snv_dict80, 0.95)

    variant_bin_and_write_to_file(snv_dict95, 'v95', gvf_set, low_conf_kmer_set, ref_file_path, out_directory_path)
    variant_bin_and_write_to_file(snv_dict80, 'v80', gvf_set, low_conf_kmer_set, ref_file_path, out_directory_path)
    variant_bin_and_write_to_file(snv_dict70, 'v70', gvf_set, low_conf_kmer_set, ref_file_path, out_directory_path)
    variant_bin_and_write_to_file(snv_dict60, 'v60', gvf_set, low_conf_kmer_set, ref_file_path, out_directory_path)
    variant_bin_and_write_to_file(snv_dict50, 'v50', gvf_set, low_conf_kmer_set, ref_file_path, out_directory_path)
    variant_bin_and_write_to_file(snv_dict40, 'v40', gvf_set, low_conf_kmer_set, ref_file_path, out_directory_path)
    variant_bin_and_write_to_file(snv_dict30, 'v30', gvf_set, low_conf_kmer_set, ref_file_path, out_directory_path)

    
if __name__ == "__main__":
    t1 = time.time()
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    t2 = time.time()
    print(f"Total time elapse {round((t2-t1)/3600,2)} hrs.")
