#!/usr/bin/python3

### Made by xiaxingquan
### April 2024

# coding = utf-8

'''
This file is some intermediate results or functions for file preprocessing
'''

import pandas as pd
import numpy as np

here = '/data/xiaxq/topic_PM1/topic_PM1_code/database/'
def add_Sumbiter():
    f = open(here + 'tmp/hot_cold_spot_modify(no-mutation-ratio)(2-sigama).txt',
             'r')
    df = pd.read_csv(here + 'tmp/final_variant-add-NumberSubmitters.txt',
                     sep='\t', low_memory=False)
    df_BLB = pd.read_csv(here + 'tmp/final_BLB_variant-add-NumberSubmitters.txt',
                         sep='\t', low_memory=False)

    r = open(here + 'tmp/test-Number-submitter.txt', 'w')

    count = 0
    for line in f:
        if count == 0:
            r.write(line.replace('\n', '') + '\t' + 'm-t-t-ratio' + '\n')
        else:
            data = line.split('\t')
            p = {}
            if data[0] == 'hotspot':
                goal = df[df['genename'] == data[2]]
                for i in range(len(goal)):
                    pos = goal.iloc[i]['aa_position']

                    number = goal.iloc[i]['NumberSubmitters']
                    if int(pos) >= int(data[4]) and int(pos) <= int(data[5]):
                        if int(pos) in p:
                            p[int(pos)] += int(number)
                        else:
                            p[int(pos)] = int(number)
                hspot = 0
                for elem in p:
                    if p[elem] >= 3:
                        hspot += 1
                l = len(p)
                if l == 0:
                    l = 1
                ratio = hspot / l
                r.write(line.replace('\n', '') + '\t' + str(ratio) + '\n')
            else:
                goal = df_BLB[df_BLB['genename'] == data[2]]
                for i in range(len(goal)):
                    pos = goal.iloc[i]['aa_position']

                    number = goal.iloc[i]['NumberSubmitters']
                    if int(pos) >= int(data[4]) and int(pos) <= int(data[5]):
                        if int(pos) in p:
                            p[int(pos)] += int(number)
                        else:
                            p[int(pos)] = int(number)
                hspot = 0
                for elem in p:
                    if p[elem] >= 3:
                        hspot += 1
                l = len(p)
                if l == 0:
                    l = 1
                ratio = hspot / l
                r.write(line.replace('\n', '') + '\t' + str(ratio) + '\n')
        count += 1

    print('ok')

    all_vcf = pd.read_csv(here + 'clinvar_20240407_temp.vcf', sep='\t',
                          low_memory=False)
    f = open(here + 'tmp/interval_first_aapos.txt', 'r')

    result = open(here + 'tmp/final_variant_add_rsid.txt', 'w')

    # PLP     1       1041335 1041346 AGCTCCTGCGCC    A       AGRN    NM_198576.4     exon5   298     hotspot 262     2224
    result.write('#CHROM' + '\t' + 'POS' + '\t' + 'ID'
                 + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'QUAL' + '\t'
                 + 'FILTER' + '\t' + 'INFO' + '\n')

    count = 0
    for line in f:
        data = line.split('\t')
        chrom = data[1]
        if chrom == 'X':
            chrom = 23
        elif chrom == 'Y':
            chrom = 24
        goal = np.array(all_vcf[all_vcf['POS'] == int(data[2])])
        for elem in goal:
            # print(elem[2])
            if elem[0] == 'X':
                elem[0] = 23
            elif elem[0] == 'Y':
                elem[0] = 24

            if int(elem[0]) == int(chrom) and elem[4] == data[5] and elem[3] == data[4]:
                ID = 'rs' + str(elem[2])
                result.write(str(data[1]) + '\t' + str(data[2]) + '\t' + ID
                             + '\t' + str(data[4]) + '\t' + str(data[5]) + '\t' + '.' + '\t'
                             + '.' + '\t' + str(data[6]) + '|' + str(data[7]) + '|'
                             + str(data[8]) + '|' + str(data[9]) + '|' + str(data[10])
                             + '|' + str(data[11]) + '|' + str(data[12]))
                break
        count += 1

    f.close()

# 1	1041335	rs3010684	AGCTCCTGCGCC	A	.	.	AGRN|NM_198576.4|exon5|298|hotspot|262|2224;CSQ=-|frameshift_variant|Transcript|ENST00000379370|5/36,-|non_coding_transcript_exon_variant|Transcript|ENST00000469403|3/3,-|upstream_gene_variant|Transcript|ENST00000479707|,-|frameshift_variant|Transcript|ENST00000620552|5/39,-|frameshift_variant|Transcript|ENST00000651234|4/38,-|frameshift_variant|Transcript|ENST00000652369|4/35,-|frameshift_variant|Transcript|NM_001305275.2|5/39,-|frameshift_variant|Transcript|NM_001364727.2|4/36,-|frameshift_variant|Transcript|NM_198576.4|5/36,-|frameshift_variant|Transcript|XM_005244749.4|5/37,-|frameshift_variant|Transcript|XM_011541429.3|5/37,-|frameshift_variant|Transcript|XM_047419836.1|5/36,-|upstream_gene_variant|Transcript|XM_047419837.1|,-|upstream_gene_variant|Transcript|XM_047419838.1|


def ENST_to_NCBI():

    human =pd.read_csv('/data/xiaxq/anaconda3/envs/vep_110/lib/ensembl-vep/database/ENST_to_NCBI.txt', sep='\t', low_memory=False)
    vep_result = open('/data/xiaxq/anaconda3/envs/vep_110/lib/ensembl-vep/result/test_240712.txt','r')

    result = open('/data/xiaxq/anaconda3/envs/vep_110/lib/ensembl-vep/result/fin_result.txt','w')

    result.write('type' + '\t' + 'chr' + '\t' + 'gene' + '\t' + 'iso' +'\t' + 'exon' +'\t' + 'aa_pos' + '\t' + 'aa_start' + '\t' + 'aa_end' + '\t' + 'exon_pos' + '\n' )

    for line in vep_result:
        if '#' not in line:
            data = line.split('\t')
            NCBI_ISO = data[7].split(';')[0].split('|')[1].split('.')[0]
            Ensembl_ISO = np.array(human[human['RNA_nucleotide_accession.version'].str.contains(NCBI_ISO)]['Ensembl_rna_identifier'])
            if len(Ensembl_ISO) == 0:
                result.write(data[7].split(';')[0].split('|')[4] + '\t' + str(data[0])  + '\t' + data[7].split(';')[0].split('|')[0] + '\t' + data[7].split(';')[0].split('|')[1] + '\t' + data[7].split(';')[0].split('|')[2] + '\t' +  data[7].split(';')[0].split('|')[3] + '\t' + data[7].split(';')[0].split('|')[5] + '\t' + data[7].split(';')[0].split('|')[6] + '\t' + 'NAN' + '\n')
                continue
            ISO_INFO = data[7].split(';')[1].split('=')[1].split(',')
            flag = 1
            for info in ISO_INFO:
                if info.split('|')[3] in Ensembl_ISO[0]:
                    flag = 0
                    result.write(data[7].split(';')[0].split('|')[4] + '\t' + str(data[0])  + '\t' + data[7].split(';')[0].split('|')[0] + '\t' + data[7].split(';')[0].split('|')[1] + '\t' + data[7].split(';')[0].split('|')[2] + '\t' +  data[7].split(';')[0].split('|')[3] + '\t' + data[7].split(';')[0].split('|')[5] + '\t' + data[7].split(';')[0].split('|')[6] + '\t' + info.split('|')[4] + '\n')

            if flag == 1:
                result.write(data[7].split(';')[0].split('|')[4] + '\t' + str(data[0])  + '\t' + data[7].split(';')[0].split('|')[0] + '\t' + data[7].split(';')[0].split('|')[1] + '\t' + data[7].split(';')[0].split('|')[2] + '\t' +  data[7].split(';')[0].split('|')[3] + '\t' + data[7].split(';')[0].split('|')[5] + '\t' + data[7].split(';')[0].split('|')[6] + '\t' + 'NAN' + '\n')

    vep_result.close()
    result.close()
    print('ok')

def var_in_hotspot_or_coldspot():

    df = pd.read_csv(
        here + 'update_database/hot_cold_result-and-profile-coefficient(1-sigama).txt',
        sep='\t', low_memory=False)
    fin_vus = open(here + 'tmp/final_variant_20240711.txt', 'r')
    coldspot = df[df['type'] == 'coldspot']
    # result = open(here + 'update_database/all_clinvar_var_in_hotspot.txt','w')
    count = 0

    for line in fin_vus:

        data = line.split('\t')
        gene = data[6]
        aa_pos = data[9]
        get_interval = np.array(coldspot[coldspot['gene'] == gene])
        if len(get_interval) > 0:
            for elem in get_interval:
                if int(aa_pos) >= int(elem[4]) and int(aa_pos) <= int(elem[5]):
                    # result.write(line)
                    count += 1
        # count += 1

    print(count)

import tabix
def compute_var_Alphamiessense():

    filename = here + 'Alphamiessense/AlphaMissense_hg38_sort.vcf.gz'
    tb = tabix.open(filename)


    vus = open(here + 'tmp/final_BLB_variant_0711.txt', 'r')
    count = 0
    result = open(here + 'BLB_Alphamiessense.txt', 'w')
    for line in vus:
        data = line.split('\t')
        records = tb.query(str(data[1]), int(data[2]) - 1, int(data[2]))
        for record in records:
            if record[3] == data[5]:
                result.write(line.replace('\n','\t') + str(record[4]) + '\n' )
        count += 1

    print('ok')

def compute_var_Mutscore():

    filename = here + 'hg38_mutscore_sort.txt.gz'
    tb = tabix.open(filename)

    vus = open(here + 'temp/var_in_hotspot.txt', 'r')
    count = 0
    result = open(here + 'temp/vus_hot_mutscore_ADD_PLP.txt', 'w')
    for line in vus:
        data = line.split('\t')
        records = tb.query(str(data[1]), int(data[2])-1, int(data[3]))
        for record in records:
            if record[4] == data[5]:
                result.write(str(record[0]) + '\t' + str(record[1]) + '\t' + str(record[2]) + '\t' + str(record[3]) + '\t' + str(record[4]) + '\t' + str(record[5]) + '\n' )
        count += 1

    print('ok')



'''
get all aa length
'''
def get_all_BLB_or_PLP_gene():
    PLP = open(here + 'update_database/final_BLB_variant.txt','r')
    file = open(here + 'update_database/all_BLB_gene.vcf','w')

    gene = []
    for line in PLP:
        genename = line.split('\t')[6]
        if genename not in gene:
            gene.append(genename)
            file.write(line)

    file.close()


def var_in_hotspot_add_rsid():
    all_vcf =pd.read_csv(here + 'clinvar_20240407.vcf', sep='\t', low_memory=False)
    f = open(here + 'update_database/all_clinvar_var_in_hotspot.txt','r')

    result = open(here + 'update_database/all_clinvar_var_in_hotspot_add_rsid.txt','w')

    #PLP     1       1041335 1041346 AGCTCCTGCGCC    A       AGRN    NM_198576.4     exon5   298
    result.write('#CHROM' + '\t' +  'POS' + '\t' + 'ID'
                 + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'QUAL' + '\t'
                 + 'FILTER' + '\t' + 'INFO' + '\n')

    count = 0
    for line in f:
        data = line.split('\t')
        chrom = data[1]
        if chrom == 'X':
            chrom = 23
        elif chrom == 'Y':
            chrom = 24
        goal = np.array(all_vcf[all_vcf['POS'] == int(data[2])])
        for elem in goal:
            #print(elem[2])
            if elem[0] == 'X':
                elem[0] = 23
            elif elem[0] == 'Y':
                elem[0] = 24


            if int(elem[0]) == int(chrom) and elem[4] == data[5] and elem[3] == data[4]:
                ID = 'rs' + str(elem[2])
                result.write(str(data[1]) + '\t' +  str(data[2]) + '\t' + ID
                 + '\t' + str(data[4]) + '\t' + str(data[5]) + '\t' + '.' + '\t'
                 + '.' + '\t' + str(data[6]) + '|' + str(data[7]) + '|'
                 + str(data[8]) + '|' + str(data[9]))
                break
        count += 1

    f.close()


#### need vep

def vep_result_process():
    f = open("/data/xiaxq/anaconda3/envs/vep_110/lib/ensembl-vep/result/test_BLB_240805.txt",'r')
    result = open("/data/xiaxq/anaconda3/envs/vep_110/lib/ensembl-vep/result/test_BLB_240805_process.txt",'w')

    result.write('#CHROM' + '\t' +  'POS' + '\t' + 'ID'
                 + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'QUAL' + '\t'
                 + 'FILTER' + '\t' + 'INFO' + '\n')
    for line in f:
        if line[0] != '#':
            data = line.split('\t')
            CSQ = data[7].split(';')[1].split('=')[1].split(',')
            iso = data[7].split(';')[0].split('|')[1]
            flag = 0
            for elem in CSQ:
                if iso.split('.')[0] in elem:
                    flag = 1
                    result.write(line.split(';')[0] + ';' + elem.replace('\n','') + '\n')
            if flag == 0:
                print(data[7].split(';')[0])

    f.close()
    result.close()

    print('ok')

def vep_result_process_2():
    f = open("/data/xiaxq/anaconda3/envs/vep_110/lib/ensembl-vep/result/test_BLB_240805_process.txt",'r')
    result = open("/data/xiaxq/anaconda3/envs/vep_110/lib/ensembl-vep/result/test_BLB_240805_process_process.txt",'w')

    result.write('#CHROM' + '\t' +  'POS' + '\t' + 'ID'
                 + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'QUAL' + '\t'
                 + 'FILTER' + '\t' + 'INFO' + '\n')
    for line in f:
        if line[0] != '#':
            data = line.split('\t')
            legth = data[7].split(';')[1].split('||||')[1].split('|')#[2].split('/')[1]
            if len(legth) > 2:
                if len(legth[2].split('/')) == 2:
                    result.write(line.split(';')[0] + '|' + legth[2].split('/')[1] + '\n')
                else:
                    print(data[7].split(';')[0])
            else:
                print(data[7].split(';')[0])

    f.close()
    result.close()

    print('ok')


def get_ALL_VUS_CON_less_2star():


    fin_vus = open(here + 'clinvar_20240407.vcf','r')

    count = 0

    result = open(here + 'update_database/ALL_VUS_CON_less2star.txt','w')

    for line in fin_vus:
        if count == 0:
            count += 1
            continue

        if 'CLNSIG=Pathogenic;' in line or 'CLNSIG=Likely_pathogenic;' in line or 'CLNSIG=Pathogenic/Likely_pathogenic;' in line or 'CLNSIG=Likely_benign;' in line or 'CLNSIG=Benign;' in line or 'CLNSIG=Benign/Likely_benign;' in line:
            if 'CLNREVSTAT=criteria_provided,_multiple_submitters' in line or 'reviewed_by_expert_panel' in line or 'CLNREVSTAT=practice_guideline' in line:
                continue
            else:
                result.write(str(line.split('\t')[0]) + '\t' + str(line.split('\t')[1]) + '\t' +  str(int(line.split('\t')[1]) + len(line.split('\t')[3]) -1) + '\t' + str(line.split('\t')[3]) + '\t' + str(line.split('\t')[4]) + '\n' )

        else:
            result.write(str(line.split('\t')[0]) + '\t' + str(line.split('\t')[1]) + '\t' +  str(int(line.split('\t')[1]) + len(line.split('\t')[3]) -1) + '\t' + str(line.split('\t')[3]) + '\t' + str(line.split('\t')[4]) + '\n' )
        count += 1

    print('ok')
