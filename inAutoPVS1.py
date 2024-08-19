import numpy as np

here = 'C:/Users/xiaxq/Desktop/topic-PM1-tep-file/'

def vaild():
    allVar_fliter = open(here + 'AutoPVS1_test.txt', 'r')
    geneList = []
    for line in allVar_fliter:
        gene = line.split('\t')[1]
        if gene not in geneList:
            geneList.append(gene)
    print(len(geneList))


def extractGene(input):
    allPLP = open(input,'r')
    geneList = []
    for line in allPLP:
        gene = line.split('\t')[2]
        if gene not in geneList:
            geneList.append(gene)

    PP2 = open(here + 'PP2.genes.hg38','r')

    PP2_PM1 = []

    for line in PP2:
        if line.replace('\n','') in geneList:
            PP2_PM1.append(line.replace('\n',''))

    print(len(PP2_PM1))


    # allVar = open(here + 'clinvar_all_2star.txt','r')
    # allVar_fliter_by_hotGene = open(here + 'clinvar_all_filter_by_hotGene_2star.txt','w')
    #
    # for line in allVar:
    #     # if line.split('\t')[-1].replace('\n', '') not in g:
    #     #     g.append(line.split('\t')[-1].replace('\n', ''))
    #     if line.split('\t')[-1].replace('\n','') in geneList:
    #         allVar_fliter_by_hotGene.write(line)
    # print('ok')

def contained_in_bed(bed_dict, chrom, start, end):
    """
    Determine whether a variant is contained in the bed region.
    :param bed_dict:
    :param chrom:
    :param start:
    :param end:
    :return: Boolean value
    """
    chrom = str(chrom)
    if "chr" not in chrom:
        chrom = "chr" + chrom
    # max(start1, start2) < min(end1, end2)
    for key in bed_dict:
        if bed_dict[key]["chrom"] == chrom and \
                max(bed_dict[key]["start"], start) < min(bed_dict[key]["end"], end):
            return True, key
    return False

def file_process(file_input):
    input = open(file_input,'r')
    hot = open(here + 'expert_curated_domains_hg38.bed','r')
    count = 0
    p = 0
    b = 0
    v = 0
    c = 0
    for line in input:
        type = line.split('\t')[0]
        chr = line.split('\t')[1]
        start = line.split('\t')[2]
        end = line.split('\t')[3]
        chrom = str(chr)
        if "chr" not in chrom:
            chrom = "chr" + chrom
        hot.seek(0)
        for h in hot:
            if h.split('\t')[0] == chrom and max(h.split('\t')[1], start) <= min(h.split('\t')[2], end):
                count += 1
                if type == 'PLP':
                    p += 1
                if type == 'BLB':
                    b += 1
                if type == 'VUS':
                    v += 1
                if type == 'CON':
                    c += 1
                #print(True)
                break
    print(count,p,b,v,c)

if __name__ == '__main__':
    """
    clinvar_all.txt clinvar_all_2star.txt
    clinvar_all_filter.txt clinvar_all_filter_2star.txt
    clinvar_all_filter_by_hotGene.txt clinvar_all_filter_by_hotGene_2star.txt
    AutuPVS1_test.txt (no type)  var_in_hotspot.txt
    """
    #file_process(here + 'clinvar_all_filter_by_hotGene_2star.txt')
    extractGene(here + 'hot_cold_spot_(1-sigama).txt')
    #vaild()

