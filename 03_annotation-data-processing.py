#!/usr/bin/python3

# Made by xiaxingquan
# April 2024

# coding = utf-8

here = '/data/xiaxq/topic_PM1/topic_PM1_code/database/'

# get amino acid position


def get_aa_position(aachange):
    try:
        change = aachange.split('.')[1]
        result = change[1:-1]
    except Exception as e:
        print(e)

    return result


def annotation_data_processing():
    '''
    inputfile:
        refGene.exonic_variant.function file by annovar
    outputfile:
        inputfile of get-mane-iso()
    '''

    input = ['PLP', 'BLB', 'VUS', 'CON']
    # input = ['BLB']

    for var_type in input:
        input_file = here + 'clinvar_annotation/update_database-annotation-' + \
            var_type + '.refGene.exonic_variant_function'
        output_file = here + 'update_database/' + var_type + '_exonic.txt'

        file = open(input_file, 'r')
        file_result = open(output_file, 'w')

        # head
        headers = [
            'chr',
            'start',
            'end',
            'ref',
            'alt',
            'genename',
            'iso',
            'exonic',
            'aa_position']
        file_result.write(f"{'\t'.join(headers)}\n")
        count = 0
        for line in file:
            count += 1
            data = line.split('\t')
            container = [data[-5], data[-4], data[-3], data[-2], data[-1]]
            txt = f"{'\t'.join(container)}"
            txt = txt.replace('\n', '')
            if data[2][-1] != ',':
                data[2] += ','
            temp = data[2].split(',')
            '''
            if len(temp) == 1:
                right = len(temp)
            else:
                right = len(temp) - 1
            '''
            for i in range(0, len(temp) - 1):
                info = temp[i].split(':')
                if len(info) == 5:
                    container = [
                        txt,
                        info[0],
                        info[1],
                        info[2],
                        get_aa_position(
                            info[4])]
                    file_result.write(f"{'\t'.join(container)}\n")

        file.close()
        file_result.close()

    print('ok')


def get_mane_iso():
    '''
    input file :
        outputfile of annotation-data-processing,mane file
    '''

    input = ['PLP', 'BLB', 'VUS', 'CON']

    for vcf in input:
        MANE = here + 'MANE/MANE.GRCh38.v1.3.summary.txt'
        ALL_ISO_FILE = here + 'update_database/' + vcf + '_exonic.txt'

        result_file = here + 'update_database/final_' + vcf + '_variant.txt'
        file = open(result_file, 'w')

        mane = pd.read_csv(MANE, sep='\t', low_memory=False)
        all_iso = pd.read_csv(ALL_ISO_FILE, sep='\t', low_memory=False)

        # get genelist
        genename = all_iso.iloc[:, 5]
        genelist = genename.unique()
        isolist = []

        for gene in genelist:
            # print(gene)
            primary = np.array(
                all_iso[all_iso['genename'] == gene]['iso'].unique())
            select = mane[mane['symbol'] == gene]
            if len(np.array(select['RefSeq_nuc'])) == 0:
                isolist.append(
                    np.array(all_iso[all_iso['genename'] == gene]['iso'])[0])
            elif np.array(select['RefSeq_nuc'])[0].split('.')[0] not in primary:
                isolist.append(primary[0])
            else:
                isolist.append(np.array(select['RefSeq_nuc'])[0])

        for index, row in all_iso.iterrows():
            if isolist[genelist.tolist().index(row['genename'])].split('.')[
                    0] == row['iso']:
                row['iso'] = isolist[genelist.tolist().index(row['genename'])]
                for i in range(len(row) - 1):
                    file.write(str(row[i]) + '\t')
                file.write(str(row[len(row) - 1]) + '\n')

        file.close()

    # add flag in prefix of file
    # for i in input:
    #   file = open(here + 'tmp/final_hg19_' + vcf + '_variant.txt','r')
    #   result = open(here + 'tmp/final_hg19_' + vcf + '_variant.vcf','w')
    #   for line in file:
    #      result.write(i + '\t' + line)
    #  file.close()
    #  result.close()

    print('ok')


if __name__ == '__main__':
    annotation_data_processing()
    get_mane_iso()
