#!/usr/bin/python3

### Made by xiaxingquan
### May 2024

# coding = utf-8

'''
 the purpose of this code is to analyze the generated results
'''
import pandas as pd
import numpy as np

here = '/data/xiaxq/topic_PM1/topic_PM1_code/'

## Analyze the distribution of cold and hot spots

df = pd.read_csv(here + 'database/tmp/hot_cold_spot-rario-0.90-variant-0.5.txt', sep='\t', low_memory=False)

only_PLP_file = open(here + 'result-analysis/only_PLP.txt','w')
only_BLB_file = open(here + 'result-analysis/only_BLB.txt','w')
all_only_PLP_file = open(here + 'result-analysis/all_only_PLP.txt','w')
all_only_BLB_file = open(here + 'result-analysis/all_only_BLB.txt','w')

genename =  df.iloc[:, 1]

### gene number
print(len(np.array(genename.unique())))

## genes with more than 2 hot and cold spots
h_c_morethan_2 = 0
## hotspot number
hot_number = 0
only_hot = []
## coldspot number
cold_number = 0
only_cold = []
## hot and cold spot number
hot_cold_number = 0

for gene in genename.unique():
    flag = df[df['gene'] == gene]
    hot = np.array(flag[flag['type'] == 'hotspot']['type'])
    cold = np.array(flag[flag['type'] == 'coldspot']['type'])
    ## genes with more than 2 hot and cold spots
    if len(hot) >= 2 and len(cold) >= 2:
        h_c_morethan_2 += 1
    ## only hotspot
    if len(hot) > 0 and len(cold) == 0:
        hot_number += 1
        only_hot.append(gene)
    ## there are both cold and hot spots
    if len(hot) > 0 and len(cold) > 0:
        hot_cold_number += 1
    ## only coldspot
    if len(hot) == 0 and len(cold) > 0:
        cold_number += 1
        only_cold.append(gene)
        
print('more than 2 hot and cold spots:',h_c_morethan_2)
print('hot:',hot_number)
print('cold:',cold_number)
print('hot_cold:',hot_cold_number)

## Analyze genes that are intolerant to mutations and genes that are tolerant to mutations

## get all variants
all_variant = pd.read_csv(here + 'database/tmp/final_variant.txt', sep='\t', low_memory=False)

## analyze genes that are all hot spots
analyze_hot = 0
all_PLP = 0
for gene in only_hot:
    goal_gene = all_variant[all_variant['genename'] == gene]
    get_PLP = goal_gene[goal_gene['type'] == 'PLP']
    get_BLB = goal_gene[goal_gene['type'] == 'BLB']
    if len(np.array(get_BLB))/(len(np.array(get_BLB)) + len(np.array(get_PLP))) <= 0.1:
        analyze_hot += 1
        only_PLP_file.write(gene)
    if len(np.array(get_BLB)) == 0:
        all_PLP += 1
        only_PLP_file.write('\t' + '100%')
    only_PLP_file.write('\n')
    #print('gene:',gene,'PLP_number:',len(np.array(get_PLP)),'BLB_number:',len(np.array(get_BLB)))

print('Genes with a PLP ratio greater than or equal to 90%:',analyze_hot)
print('Genes with only PLP:',all_PLP)

## analyze genes that are all cold spots
analyze_cold = 0
all_BLB = 0
for gene in only_cold:
    goal_gene = all_variant[all_variant['genename'] == gene]
    get_PLP = goal_gene[goal_gene['type'] == 'PLP']
    get_BLB = goal_gene[goal_gene['type'] == 'BLB']
    if len(np.array(get_PLP))/(len(np.array(get_BLB)) + len(np.array(get_PLP))) <= 0.1:
        analyze_cold += 1
        only_BLB_file.write(gene)
    if len(np.array(get_PLP)) == 0:
        all_BLB += 1
        only_BLB_file.write('\t' + '100%')
    #print('gene:',gene,'PLP_number:',len(np.array(get_PLP)),'BLB_number:',len(np.array(get_BLB)))
    only_BLB_file.write('\n')
   
print('Genes with a BLB ratio greater than or equal to 90%:',analyze_cold)
print('Genes with only BLB:',all_BLB)
## all gene
all_gene = all_variant['genename'].unique()
analyze_cold = 0
analyze_hot = 0
all_PLP = 0
all_BLB = 0

print('number of all gene:',len(all_gene))

for gene in all_gene:
    goal_gene = all_variant[all_variant['genename'] == gene]
    get_PLP = goal_gene[goal_gene['type'] == 'PLP']
    get_BLB = goal_gene[goal_gene['type'] == 'BLB']
    if len(np.array(get_BLB)) + len(np.array(get_PLP)) > 0:
        if len(np.array(get_PLP))/(len(np.array(get_BLB)) + len(np.array(get_PLP))) <= 0.1:
            analyze_cold += 1
            all_only_BLB_file.write(gene)
        if len(np.array(get_PLP)) == 0:
            all_BLB += 1
            all_only_BLB_file.write('\t' + '100%')
        #print('gene:',gene,'PLP_number:',len(np.array(get_PLP)),'BLB_number:',len(np.array(get_BLB)))
        if len(np.array(get_BLB))/(len(np.array(get_BLB)) + len(np.array(get_PLP))) <= 0.1:
            analyze_hot += 1
            all_only_PLP_file.write(gene)
        if len(np.array(get_BLB)) == 0:
            all_BLB += 1
            all_only_PLP_file.write('\t' + '100%')
     
only_PLP_file.close()
only_BLB_file.close()
all_only_PLP_file.close()
all_only_BLB_file.close()
       
print('ALL genes with a BLB ratio greater than or equal to 90%:',analyze_cold)
print('ALL genes with a PLP ratio greater than or equal to 90%:',analyze_hot)
print('ALL genes with only BLB:',all_BLB)
print('ALL genes with only PLP:',all_PLP)