#!/usr/bin/python3

# Made by xiaxingquan
# April 2024

# coding = utf-8

"""
    input file :
        outputfile of EM.py
    output file:
        final result of hotspot area and coldspot area

    The purpose of this program is to process the initial results,
    screen areas with a ratio greater than or equal to 0.90, and obtain the final result
"""

here = '/data/xiaxq/topic_PM1/topic_PM1_code/database/'


def initial_result_processing():

    file = open(here + 'tmp/Initial_results_modify(2-sigama)', 'r')
    result = open(
        here + 'tmp/hot_cold_spot_modify(no-mutation-ratio)(2-sigama).txt', 'w')

    headers = [
        'type',
        'chr',
        'gene',
        'iso',
        'aa_start_pos',
        'aa_end_pos',
        'PLP_number',
        'BLB_number',
        'Ratio_ratio']
    result.write(f"{'\t'.join(headers)}\n")

    for line in file:
        try:
            data = line.split('\t')
            chromosome = data[0]
            gene = data[1]
            iso = data[2]
            # hot spot area
            PLP = data[3].replace('(PLP)', '')
            if len(PLP) != 0:
                hotspot = PLP.split(';')
                for i in range(0, len(hotspot) - 1):
                    content = hotspot[i].split(':')
                    if content[0] == '0':
                        content[0] = '1'
                    # if float(content[4]) >= 0.90 and float(content[2]) /
                    # (float(content[1]) - float(content[0])) >= 0.5:
                    if float(content[4]) >= 0.90:
                        container = ['hotspot', chromosome, gene, iso,
                                     content[0], content[1], content[2], content[3],
                                     content[4]]
                        result.write(f"{'\t'.join(container)}\n")

            # cold spot area
            BLB = data[4].replace('(BLB)', '')
            if len(BLB) != 0:
                coldspot = BLB.split(';')
                for i in range(0, len(coldspot) - 1):
                    content = coldspot[i].split(':')
                    if content[0] == '0':
                        content[0] = '1'
                    # if float(content[4]) >= 0.90 and float(content[3]) /
                    # (float(content[1]) - float(content[0])) >= 0.5:
                    if float(content[4]) >= 0.90:
                        container = ['coldspot', chromosome, gene, iso,
                                     content[0], content[1], content[2], content[3],
                                     content[4]]
                        result.write(f"{'\t'.join(container)}\n")
        except Exception as e:
            print(e)

    print('ok')
    file.close()
    result.close()


if __name__ == '__main__':
    initial_result_processing()
