#!/usr/bin/python3

### Made by xiaxingquan
### May 2024

# coding = utf-8


###The purpose of this program is to automatically score user input mutations. 
###A score of 0 indicates that the mutation is not in the clustering area, 
### 1 indicates that the mutation is in the hot spot area, and -1 indicates that the mutation is in the cold spot area.


import subprocess
import argparse


here = '/data/xiaxq/topic_PM1/topic_PM1_code/'

parser = argparse.ArgumentParser()
## input file
parser.add_argument('-i', type=str, default = None)
## output file
parser.add_argument('-o', type=str, default = None)
## human genome version
parser.add_argument('-buildver', type=str, default= 'hg38')
args = parser.parse_args()


## hg19 version to hg38
def hg19_to_hg38(ifile,ofile,tranfomertype):
    #Define the command-line parameters for tranformer
    tranformer_command = [
        'sh',
        here + 'hg19-to-hg38.sh',
        ifile,
        ofile,
        tranfomertype
    ]
    
    #run command
    subprocess.run(tranformer_command)

## process the path of input files
if args.i[0] != '/' and args.i[0] != '.':
    inputfile = here + args.i
else:
    inputfile = args.i
    
## process hg19 or hg38
if args.buildver == 'hg19':
    final_inputfile = here + 'database/tmp/hg19-to-hg38-temp.txt'
    hg19_to_hg38(inputfile,final_inputfile,'hg19Tohg38')
else:
    final_inputfile = inputfile
    

## process the path of output files
if args.o[0] != '/' and args.o[0] != '.':
    outfile = here + args.o
else:
    outfile = args.o

#input_vcf, output_vcf, reference_genome, annotation_database

print('NOTICE: Annotating input file by annovar...')

def run_annovar():
    #Define the command-line parameters for annovar
    ANNOVAR_command = [
        'perl',
        here + 'annovar/annotate_variation.pl', 
        '-geneanno', 
        '-dbtype', 'refGene',
        '-outfile', here + 'annovar/result/Auto-annotation-annovar',
        '-buildver', 'hg38',
        '-exonsort', '-nofirstcodondel',
        final_inputfile, 
        here + 'annovar/humandb/'
    ]
    
    #run annovar
    subprocess.run(ANNOVAR_command)


'''input_vcf = 'input.vcf'
output_vcf = 'output.vcf'
reference_genome = 'GRCh38'
annotation_database = 'dbNSFP,gnomAD'
plugin_parameter = 'my_plugin_param_value'''


run_annovar()

print('NOTICE: Annotating by annovar result write to ' + here + 'annovar/result/Auto-annotation-annovar.exonic_variant_function.txt')

import pandas as pd
import numpy as np


###get amino acid position
def get_aa_position(aachange):
    result = ''
    flag = False
    for i in aachange:
        if i in '0123456789':
            flag = True
            result = result + i
        else:
            if flag:
                break

    return result

### all hot-cold spot file
hot_cold_result_path  = here + 'database/tmp/hot_cold_spot-rario-0.90-variant-0.5.txt'

### input files that users need to annotate
input_file_path = here + 'annovar/result/Auto-annotation-annovar.exonic_variant_function'

### output file of result
result_file = outfile

print('NOTICE: Annotating input file...')

file = open(result_file,'w')

### head
file.write('chr' + '\t' + 'start' + '\t' + 'end' + '\t' + 'alt' + '\t' + 'ref' + '\t' + 'otherInfo' + '\t' + 'hot_cold_score' + '\n')

hot_cold_result = pd.read_csv(hot_cold_result_path, sep='\t', low_memory=False)
input_file = open(input_file_path,'r')

### traverse all variants
for line in input_file:
    # the initialization score is 0, which means that the variant is initialized as neither a cold spot nor a hot spot
    score = 0
    data = line.split('\t')
    # store variant location information
    pos_info = data[2].split(',')
    # variant genename
    gene = pos_info[0].split(':')[0]
    # find the gene where the current variant is located in all cold hotspots
    annotation_gene = hot_cold_result[hot_cold_result['gene'] == gene]
    # find isoform
    annotation_iso = np.array(annotation_gene['iso'])
    # determine whether the gene has a cold or hot spot
    if len(annotation_iso) != 0:
        iso = annotation_iso[0]
        # record interval
        annotation_info = []
        for elem in np.array(annotation_gene):
            # if coldspot area
            if elem[0] == 'coldspot':
                annotation_info.append([-1,elem[3],elem[4]])
            else:
                annotation_info.append([1,elem[3],elem[4]])
        for elem in pos_info:
            info = elem.split(':')
            # determine whether the transcript matches
            if info[1] == annotation_iso[0].split('.')[0] and len(info) == 5:
                aa_pos = get_aa_position(info[4])
                for i in annotation_info:

                    ## determine if the variant is within the interval
                    if int(aa_pos) >= i[1] and int(aa_pos) <= i[2]:
                        score = i[0]
    # write result to result file
    file.write(data[3] + '\t' + data[4] + '\t' + data[5] + '\t' + data[6] + '\t' + data[7] + '\t' + data[8].replace('\n','') + '\t' + str(score) + '\n' )


file.close()
input_file.close()

import os

if args.buildver == 'hg19':
    
    ## transformer hg38 to hg19
    hg19_to_hg38(result_file,here + 'xiaxq','hg38Tohg19')

    ## hg19 result process
    R_file = open(result_file,'w')
    F_file = open(here + 'xiaxq','r')
    R_file.write('chr' + '\t' + 'start' + '\t' + 'end' + '\t' + 'alt' + '\t' + 'ref' + '\t' + 'otherInfo' + '\t' + 'hot_cold_score' + '\n')
    for line in F_file:
        R_file.write(line)
    
    R_file.close()
    F_file.close()
    
    os.remove(here + 'xiaxq')

    os.remove(here + 'database/tmp/hg19-to-hg38-temp.txt')

print('NOTICE: Annotating result write to ' + outfile)
print('NOTICE:Annotation completed')
