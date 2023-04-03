#!/usr/bin/env python
# coding: utf-8

import os
import sys

def main(bed_directories):

    v30_set = set()
    v40_set = set()
    v50_set = set()
    v60_set = set()
    v70_set = set()
    v80_set = set()
    v95_set = set()




    for directory in bed_directories:
        v30 = open(f'{directory}/v30_pooled_{directory.split("/")[-1].split("_")[0]}.bed', 'w')
        v40 = open(f'{directory}/v40_pooled_{directory.split("/")[-1].split("_")[0]}.bed', 'w')
        v50 = open(f'{directory}/v50_pooled_{directory.split("/")[-1].split("_")[0]}.bed', 'w')
        v60 = open(f'{directory}/v60_pooled_{directory.split("/")[-1].split("_")[0]}.bed', 'w')
        v70 = open(f'{directory}/v70_pooled_{directory.split("/")[-1].split("_")[0]}.bed', 'w')
        v80 = open(f'{directory}/v80_pooled_{directory.split("/")[-1].split("_")[0]}.bed', 'w')
        v95 = open(f'{directory}/v95_pooled_{directory.split("/")[-1].split("_")[0]}.bed', 'w')
        for file in os.listdir(directory):
            if 'pooled' in file:
                continue
            read_file = open(f'{directory}/{file}', 'r')
            for line in read_file:
                if 'v30' in file:
                    v30.write(line)
                    v30_set.add(line)
                if 'v40' in file:
                    v40.write(line)
                    v40_set.add(line)
                if 'v50' in file:
                    v50.write(line)
                    v50_set.add(line)
                if 'v60' in file:
                    v60.write(line)
                    v60_set.add(line)
                if 'v70' in file:
                    v70.write(line)
                    v70_set.add(line)
                if 'v80' in file:
                    v80.write(line)
                    v80_set.add(line)
                if 'v95' in file:
                    v95.write(line)
                    v95_set.add(line)
        v30.close()
        v40.close()
        v50.close()
        v60.close()
        v70.close()
        v80.close()
        v95.close()


    v30_pan = open('./pooled_bed/v30_pan.bed', 'w')
    for line in v30_set:
        v30_pan.write(line)
    v30_pan.close()


    v40_pan = open('./pooled_bed/v40_pan.bed', 'w')
    for line in v40_set:
        v40_pan.write(line)
    v40_pan.close()


    v50_pan = open('./pooled_bed/v50_pan.bed', 'w')
    for line in v50_set:
        v50_pan.write(line)
    v50_pan.close()


    v60_pan = open('./pooled_bed/v60_pan.bed', 'w')
    for line in v60_set:
        v60_pan.write(line)
    v60_pan.close()


    v70_pan = open('./pooled_bed/v70_pan.bed', 'w')
    for line in v70_set:
        v70_pan.write(line)
    v70_pan.close()


    v80_pan = open('./pooled_bed/v80_pan.bed', 'w')
    for line in v80_set:
        v80_pan.write(line)
    v80_pan.close()


    v95_pan = open('./pooled_bed/v95_pan.bed', 'w')
    for line in v95_set:
        v95_pan.write(line)
    v95_pan.close()

if __name__ == "__main__":
    b_direct = []
    for arg in sys.argv[1:]:
        b_direct.append(arg)
    main(b_direct)

