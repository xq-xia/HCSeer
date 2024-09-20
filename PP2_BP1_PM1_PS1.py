import numpy as np
import pandas as pd

'''
PS1 and PM1
PP2 and PM1
PM5 and PM1
coldspot and BP1
'''

here = 'C:/Users/xiaxq/Desktop/topic-PM1-tep-file/'

def process():
    PLP = open(here + 'update_database-annotation-PLP.refGene.exonic_variant_function','r')
    BLB = open(here + 'update_database-annotation-BLB.refGene.exonic_variant_function','r')
    result = open(here + 'PP2_BP1_process.txt','w')

    ##得到所有PLP基因列表
    PLP_geneList = []
    for line in PLP:
        if line.split('\t')[2].split(':')[0] not in PLP_geneList:
            PLP_geneList.append(line.split('\t')[2].split(':')[0])

    ##文件指针归0
    PLP.seek(0)

    count = 0
    resultList = []
    for gene in PLP_geneList:
        count += 1
        ###基因名、PLP无义数量、PLP移码数量、PLP错义数量、BLB无义数量、BLB移码数量、BLB错义数量
        resultList.append([gene,0,0,0,0,0,0])
        PLP.seek(0)
        for plp in PLP:
            geneName = plp.split('\t')[2].split(':')[0]
            type = plp.split('\t')[1]
            if gene == geneName:
                if 'stopgain' in type:
                    resultList[-1][1] += 1
                elif 'frameshift' in type:
                    resultList[-1][2] += 1
                else:
                    resultList[-1][3] += 1

        BLB.seek(0)
        for blb in BLB:
            geneName = blb.split('\t')[2].split(':')[0]
            type = blb.split('\t')[1]
            if gene == geneName:
                if 'stopgain' in type:
                    resultList[-1][4] += 1
                elif 'frameshift' in type:
                    resultList[-1][5] += 1
                else:
                    resultList[-1][6] += 1
        print(f'{count} genes have been completed')
    for re in resultList:
        result.write(re[0] + '\t' + str(re[1]) + '\t' + str(re[2]) + '\t' + str(re[3]) + '\t' + str(re[4]) + '\t' + str(re[5])  + '\t' + str(re[6]) + '\n')

    PLP.close()
    BLB.close()
    result.close()
    print('ok')

def PP2():
    primary = open(here + 'PP2_BP1_process.txt','r')
    PP2_gene = []

    for line in primary:
        data = line.split('\t')
        pp2_t1,pp2_t2 =0,0

        if int(data[3]) / (int(data[1]) + int(data[2]) + int(data[3])) >= 0.8:
            pp2_t1 = 1
        if int(data[-1]) /(int(data[3]) + int(data[-1]) + 0.01) <= 0.1:
            pp2_t2 = 1

        if pp2_t2 == 1 and pp2_t2 == 1:
            PP2_gene.append(data[0])

    primary.close()
    print(len(PP2_gene))
    return PP2_gene

def BP1():
    primary = open(here + 'PP2_BP1_process.txt', 'r')
    BP1_gene = []

    for line in primary:
        data = line.split('\t')
        if int(data[3]) / (int(data[1]) + int(data[2]) + int(data[3])) <= 0.1:
            BP1_gene.append(data[0])

    primary.close()
    return BP1_gene

def PM1_PP2():
    PP2_gene = PP2()
    df = pd.read_csv(here + 'hot_cold_result-and-profile-coefficient(1-sigama).txt',
                     sep='\t', low_memory=False)
    hotspot = df[df['type'] == 'hotspot']
    hotspotGene = hotspot['gene'].unique()

    for elem in PP2_gene:
        if elem in hotspotGene:
            print(elem)


def coldspot_BP1():
    BP1_gene = BP1()
    print('all BP1 gene:',len(BP1_gene))
    df = pd.read_csv(here + 'hot_cold_result-and-profile-coefficient(1-sigama).txt',
                     sep='\t', low_memory=False)
    coldspot = df[df['type'] == 'coldspot']
    coldspotGene = coldspot['gene'].unique()

    count = 0
    for elem in BP1_gene:
        if elem in coldspotGene:
            print(elem)
            count += 1

    print('coldspot and BP1:',count)


def PS1_PM5():
    all_VUS = pd.read_csv(here + 'final_ALL_VUS_CON_less2star_update_variant.txt',sep='\t',low_memory=False)
    all_PLP = open(here + 'final_PLPupdate_variant(compute-PS1).txt', 'r')

    PS1 = open(here + 'PS1_variant.txt', 'w')
    PM5 = open(here + 'PM5_variant.txt', 'w')
    #PS1 = []
    #PM5 = []

    count = 0
    for var in all_PLP:
        data = var.split('\t')
        aa_pos = int(data[-1].replace('\n',''))
        gene = data[7]
        ref_aa = data[5]
        alt_aa = data[6]

        get_gene = np.array(all_VUS[all_VUS['genename'] == gene])
        for elem in get_gene:
            if elem[-1] == aa_pos:
                if elem[5] == ref_aa and elem[6] == alt_aa:
                    PS1.write(elem[0] + '\t' + str(elem[1]) + '\t'
                              + str(elem[2]) + '\t' + str(elem[3])
                              + '\t' + str(elem[4]) + '\t'
                              + str(elem[5]) + '\t' + str(elem[6])
                              + '\t' + str(elem[7]) + '\t' + str(elem[8])
                              + '\t' + str(elem[9]) + '\t' + str(elem[10]) + '\n')

                if elem[5] == ref_aa and elem[6] != alt_aa:
                    PM5.write(elem[0] + '\t' + str(elem[1]) + '\t'
                              + str(elem[2]) + '\t' + str(elem[3])
                              + '\t' + str(elem[4]) + '\t'
                              + str(elem[5]) + '\t' + str(elem[6])
                              + '\t' + str(elem[7]) + '\t' + str(elem[8])
                              + '\t' + str(elem[9]) + '\t' + str(elem[10]) + '\n')

        count += 1

        if count % 100 == 0:
            print(f'{count} mutations have been completed')

    all_PLP.close()

    print('ok')

def PS1_PM5_filter_by_dbscsnv11():
    PS1 = open(here + 'PS1_variant.txt', 'r')
    PM5 = open(here + 'PM5_variant.txt', 'r')

    dbscnv11 = pd.read_csv(here + 'hg38_dbscsnv11.txt',
                           sep='\t', low_memory=False)

    PS1_filter = open(here + 'final_PS1_variant.txt', 'w')
    PM5_filter = open(here + 'final_PM5_variant.txt', 'w')

    threshold = 0.6

    ps1_count = 0
    pre_ps1 = []
    for ps1 in PS1:
        ps1_count += 1

        if ps1_count % 100 == 0:
            print(f'{ps1_count} PS1 mutations have been completed')
        if ps1 in pre_ps1:
            continue
        pre_ps1.append(ps1)
        data = ps1.split('\t')
        chr = data[0]
        start = data[1]
        alt = data[4]
        ref_aa = data[5]
        alt_aa = data[6]
        goal_1 = dbscnv11[dbscnv11['Start'] == int(start)]
        goal_2 = np.array(goal_1[goal_1['Alt'] == alt])

        if len(goal_2) > 0:
            for elem in goal_2:
                if str(elem[0]) == str(chr):
                    if float(elem[-1]) <= threshold and float(elem[-2]) <= threshold:
                        PS1_filter.write(ps1)
        elif ref_aa != alt_aa:
            PS1_filter.write(ps1)
        else:
            continue


    print('all PS1 have been completed')

    pm5_count = 0
    pre_pm5 = []
    for pm5 in PM5:
        pm5_count += 1

        if pm5_count % 100 == 0:
            print(f'{pm5_count} PM5 mutations have been completed')

        if pm5 in pre_pm5:
            continue
        pre_pm5.append(pm5)
        data = pm5.split('\t')
        chr = data[0]
        start = data[1]
        alt = data[4]
        ref_aa = data[5]
        alt_aa = data[6]
        goal_1 = dbscnv11[dbscnv11['Start'] == int(start)]
        goal_2 = np.array(goal_1[goal_1['Alt'] == alt])

        if len(goal_2) > 0:
            for elem in goal_2:
                if str(elem[0]) == str(chr):
                    if float(elem[-1]) <= threshold and float(elem[-2]) <= threshold:
                        PM5_filter.write(pm5)
        elif ref_aa != alt_aa:
            PM5_filter.write(pm5)
        else:
            continue

    print('all PM5 have been completed')

def PS1_in_hotspot():
    df = pd.read_csv(here + 'hot_cold_result-and-profile-coefficient(1-sigama).txt',
                     sep='\t', low_memory=False)
    PS1 = open(here + 'final_PS1_variant.txt', 'r')
    hotspot = df[df['type'] == 'hotspot']

    result = open(here + 'final_PS1_in_hotspot.txt', 'w')

    count = 0

    for line in PS1:
        data = line.split('\t')
        gene = data[7]
        aa_pos = data[-1]
        get_interval = np.array(hotspot[hotspot['gene'] == gene])
        if len(get_interval) > 0:
            for elem in get_interval:
                if int(aa_pos) >= int(elem[4]) and int(aa_pos) <= int(elem[5]):
                    result.write(line)
        count += 1

    print('ok')

def PM5_in_hotspot():
    df = pd.read_csv(here + 'hot_cold_result-and-profile-coefficient(1-sigama).txt',
                     sep='\t', low_memory=False)
    PM5 = open(here + 'final_PM5_variant.txt', 'r')
    hotspot = df[df['type'] == 'hotspot']

    result = open(here + 'final_PM5_in_hotspot.txt', 'w')

    count = 0

    for line in PM5:
        data = line.split('\t')
        gene = data[7]
        aa_pos = data[-1]
        get_interval = np.array(hotspot[hotspot['gene'] == gene])
        if len(get_interval) > 0:
            for elem in get_interval:
                if int(aa_pos) >= int(elem[4]) and int(aa_pos) <= int(elem[5]):
                    result.write(line)
        count += 1

    print('ok')

def test():
    all_VUS = pd.read_csv(here + 'final_ALL_VUS_CON_less2star_update_variant.txt',
                          sep='\t', low_memory=False)
    get_gene = np.array(all_VUS[all_VUS['genename'] == 'BRCA1'])
    for elem in get_gene:
        print(elem)

def PS1_var_in_hotspot():
    PS1 = open(here + 'PS1.AA.change.patho.hg38','r')
    clinvar_in_hotspot = open('C:/Users/xiaxq/Desktop/tmp/all_clinvar_var_in_hotspot_add_rsid.txt','r')
    PLP = open(here + 'final_PLP_variant.txt','r')
    PS1_variant_in_hotspot = open(here + 'PS1_variant_in_hotspot.txt','w')
    count = 0
    for line in PS1:
        data = line.split('\t')
        chr = data[0]
        pos = data[1]
        ref = data[3]
        alt = data[4]
        clinvar_in_hotspot.seek(0)

        for elem in clinvar_in_hotspot:
            info = elem.split('\t')
            if info[0] == chr and info[1] == pos and info[3] == ref and info[4] == alt:
                PLP.seek(0)
                flag = False
                for plp_elem in PLP:
                    plp = plp_elem.split('\t')
                    if plp[1] == chr and plp[2] == pos and plp[4] == ref and plp[5] == alt:
                        flag = True
                if flag == False:
                    print(elem.replace('\n',''))
                    PS1_variant_in_hotspot.write(elem)
                    count += 1

    print(count)
    PS1.close()
    clinvar_in_hotspot.close()

if __name__ == '__main__':
    '''
    PP2()
    PM1_PP2()
    print('===========================================')
    coldspot_BP1()
    '''
    '''
    PS1_in_hotspot()
    PM5_in_hotspot()
    '''
    coldspot_BP1()

