#!/usr/bin/python3

### Made by xiaxingquan
### june 2024

# coding = utf-8

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import fisher_exact

here = '/data/xiaxq/topic_PM1/topic_PM1_code/database/'

'''
The purpose of this code is to calculate the p-value of the hot and cold spot results on each gene
'''

def compute_p_value():

    df = pd.read_csv(
        here + 'tmp/hot_cold_spot_modify(no-mutation-ratio)(3-sigama).txt', sep='\t',
        low_memory=False)
    all_variant = pd.read_csv(here + 'tmp/final_variant.txt', sep='\t',
                              low_memory=False)
    result = open(
        here + 'p-value/cold_hotspot_p-value(no-mutation-ratio)(3-sigama).txt', 'w')
    
    result.write(
        'gene' + '\t' + 'PLP_in_hotspot' + '\t' + 'BLB_in_hotspot' + '\t' + 'PLP_in_all' + '\t' + 'BLB_in_all' + '\t' + 'P-value' + '\t' + 'oddsratio' + '\n')
    
    genename = df.iloc[:, 2]
    
    ## hotspot number
    hot_number = 0
    only_hot = []
    ## coldspot number
    cold_number = 0
    only_cold = []
    ## hot and cold spot number
    hot_cold_number = 0
    hot_and_cold = []
    
    for gene in genename.unique():
        flag = df[df['gene'] == gene]
        hot = np.array(flag[flag['type'] == 'hotspot']['type'])
        cold = np.array(flag[flag['type'] == 'coldspot']['type'])
        ## only hotspot
        if len(hot) > 0 and len(cold) == 0:
            hot_number += 1
            only_hot.append(gene)
        ## there are both cold and hot spots
        if len(hot) > 0 and len(cold) > 0:
            hot_cold_number += 1
            hot_and_cold.append(gene)
        ## only coldspot
        if len(hot) == 0 and len(cold) > 0:
            cold_number += 1
            only_cold.append(gene)
    
    print('hot:', hot_number)
    print('cold:', cold_number)
    print('hot_cold:', hot_cold_number)
    
    for gene in hot_and_cold:
        # get gene
        get_result = df[df['gene'] == gene]
        ## hotspot p-value analyse
        '''
        get_result = get_result[get_result['type'] == 'hotspot']
        #PLP number
        PLP_in_hotspot = 0
        #BLB number
        BLB_in_hotspot = 0
        for i in get_result['PLP_number']:
            PLP_in_hotspot += int(i)
        for i in get_result['BLB_number']:
            BLB_in_hotspot += int(i)
    
        #get PLP number in all variants
        PLP_in_all = 0
        #get BLB number in all variants
        BLB_in_all = 0
    
        get_gene = all_variant[all_variant['genename'] == gene]
        PLP = get_gene[get_gene['type'] == 'PLP']
        BLB = get_gene[get_gene['type'] == 'BLB']
    
        PLP_in_all = len(PLP)
        BLB_in_all = len(BLB)
    '''
    
        ## coldspot p-value analyse
    
        get_result = get_result[get_result['type'] == 'coldspot']
        # PLP number
        PLP_in_coldspot = 0
        # BLB number
        BLB_in_coldspot = 0
        for i in get_result['PLP_number']:
            PLP_in_coldspot += int(i)
        for i in get_result['BLB_number']:
            BLB_in_coldspot += int(i)
    
        # get PLP number in all variants
        PLP_in_all = 0
        # get BLB number in all variants
        BLB_in_all = 0
    
        get_gene = all_variant[all_variant['genename'] == gene]
        PLP = get_gene[get_gene['type'] == 'PLP']
        BLB = get_gene[get_gene['type'] == 'BLB']
    
        PLP_in_all = len(PLP)
        BLB_in_all = len(BLB)
    
        # Building contingency tables
        ##hotspot
        # table = [[PLP_in_hotspot, BLB_in_hotspot],  #Number of disease-related and non disease-related variants in hotspot areas
        # [PLP_in_all - PLP_in_hotspot, BLB_in_all - BLB_in_hotspot]] #Number of disease-related and non disease-related variants in non hotspot areas
    
        ##coldspot
        table = [[BLB_in_coldspot, PLP_in_coldspot],
                 # Number of disease-related and non disease-related variants in coldspot areas
                 [BLB_in_all - BLB_in_coldspot,
                  PLP_in_all - PLP_in_coldspot]]  # Number of disease-related and non disease-related variants in non coldspot areas
    
        # Perform Fisher's exact test
        oddsratio, p_value = fisher_exact(table)
        # print(gene,table)
        result.write(gene + '\t' + str(BLB_in_coldspot) + '\t' + str(PLP_in_coldspot) + '\t' + str(BLB_in_all) + '\t' + str(
            PLP_in_all) + '\t' + str(p_value) + '\t' + str(oddsratio) + '\n')
        # result.write(gene + '\t' + str(PLP_in_hotspot) + '\t' + str(BLB_in_hotspot) + '\t' + str(PLP_in_all) + '\t' + str(BLB_in_all) + '\t' + str(p_value) + '\t' + str(oddsratio) + '\n')
    
    result.close()
    print('ok')
    
if __name__ == '__main__':
    compute_p_value()
