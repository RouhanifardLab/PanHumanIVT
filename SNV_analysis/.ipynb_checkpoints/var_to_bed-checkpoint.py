import os
import sys


#Function that takes a file path, converts our output data to a bed file with the result.
def var_to_bed(file_path):
    rgb_dict = {'G' : '209,113,5', 'T' : '255,0,0', 'A' : '50,241,50', 'C' : '0,0,255'}
    write_path = file_path[:-4] + ".bed"
    read_fh = open(file_path, 'r')
    write_fh = open(write_path, 'w')
    for line in read_fh:
        if "Position" in line:
            continue
        split_line = line.split()[:3]
        chrom = split_line[0]
        start = int(split_line[1]) - 1
        stop = start + 1
        name = split_line[2]
        line = [chrom, str(start), str(stop), f"Observed Variant: {name}", "0", ".", str(start), str(stop), rgb_dict[name]]
        line = "\t".join(line)+"\n"
        write_fh.write(line)
    read_fh.close()
    write_fh.close()

def main(directory):
    files = os.listdir(directory)
    for file in files:
        if '.tsv' not in file or 'pooled' in file:
            continue
        var_to_bed(f"{directory}/{file}")

if __name__ == "__main__":
    main(sys.argv[1])




