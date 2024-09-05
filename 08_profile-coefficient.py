#!/usr/bin/python3

### Made by xiaxingquan
### june 2024

# coding = utf-8

import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
import random
from sklearn.metrics import silhouette_samples


here = '/data/xiaxq/topic_PM1/topic_PM1_code/database/'

def compute_profile_coefficient():
    all_variant = pd.read_csv(here + 'tmp/final_variant.txt', sep='\t',
                              low_memory=False)
    df = pd.read_csv(
        here + 'tmp/hot_cold_spot_modify(no-mutation-ratio)(2-sigama).txt',
        sep='\t', low_memory=False)
    result_file = open(
        here + 'tmp/profile-coefficient(no-mutation-ratio)(2-sigama).txt', 'w')

    result_file.write('gene' + '\t' + 'interval' + '\t' + 'profile-coefficient' + '\n')

    ###Calculate the contour coefficient of the hotspot interval
    hot = df[df['type'] == 'hotspot']

    ##get hotspot gene list
    hot_geneList = hot['gene'].unique()

    for gene in hot_geneList:
        ##get variant
        get_variant = all_variant[all_variant['genename'] == gene]
        PLP = get_variant[get_variant['type'] == 'PLP']

        inter = hot[hot['gene'] == gene]

        ##get interval
        intervals = []
        for i in range(0, len(list(inter['aa_start_pos']))):
            intervals.append([int(list(inter['aa_start_pos'])[i]), int(list(inter['aa_end_pos'])[i])])
            # print(list(inter['aa_start_pos']))
        samples = list(PLP['aa_position'])

        ##sort intervals
        intervals = sorted(intervals, key=lambda x: x[0])
        label = []
        new_interval = []

        count = 0
        ##Assign labels to samples based on intervals
        for i in intervals:
            if count == 0:
                if 1 < intervals[0][0]:
                    new_interval.append([1, intervals[0][0] - 1])
                    new_interval.append(intervals[0])
                else:
                    new_interval.append(intervals[0])
            else:
                new_interval.append([intervals[intervals.index(i) - 1][1] + 1, i[0] - 1])
                new_interval.append(i)

            count += 1
        for i in samples:
            for j in new_interval:
                if i >= j[0] and i <= j[1]:
                    label.append(new_interval.index(j))

            if i > new_interval[len(new_interval) - 1][1]:
                label.append(len(new_interval))

        label_interval = []
        for elem in intervals:
            label_interval.append(new_interval.index(elem))
        ##duplicate removal
        process = list(set(label))
        result = []
        for g in label_interval:
            result.append(process.index(g))

        sample = []
        for i in samples:
            sample.append([i])

        samples = np.array(sample)

        if len(list(set(label))) > 1:
            labels = np.array(label)
            sample_silhouette_values = silhouette_samples(samples, labels)

            silhouette_per_cluster = {}
            for label in np.unique(labels):
                silhouette_per_cluster[label] = np.mean(sample_silhouette_values[labels == label])

            for i in label_interval:
                result_file.write(gene + '\t' + str(new_interval[i][0]) + ',' + str(new_interval[i][1]) + '\t' + str(
                    silhouette_per_cluster[i]) + '\n')
                ## If there is only one cluster in the transcript, the contour coefficient is -1 to 1
        elif len(list(set(label))) == 1 and len(new_interval) == 1:
            result_file.write(gene + '\t' + str(new_interval[0][0]) + ',' + str(new_interval[0][1]) + '\t' + str(
                format(random.uniform(-1, 1), '.8f')) + '\n')
        else:
            result_file.write(gene + '\t' + str(new_interval[1][0]) + ',' + str(new_interval[1][1]) + '\t' + str(
                format(random.uniform(-1, 1), '.8f')) + '\n')

    ###Calculate the contour coefficient of the coldspot interval
    cold = df[df['type'] == 'coldspot']

    cold_geneList = cold['gene'].unique()

    for gene in cold_geneList:
        get_variant = all_variant[all_variant['genename'] == gene]
        BLB = get_variant[get_variant['type'] == 'BLB']

        inter = cold[cold['gene'] == gene]

        intervals = []
        for i in range(0, len(list(inter['aa_start_pos']))):
            intervals.append([int(list(inter['aa_start_pos'])[i]), int(list(inter['aa_end_pos'])[i])])
            # print(list(inter['aa_start_pos']))
        samples = list(BLB['aa_position'])

        intervals = sorted(intervals, key=lambda x: x[0])
        label = []
        new_interval = []

        count = 0
        for i in intervals:
            if count == 0:
                if 1 < intervals[0][0]:
                    new_interval.append([1, intervals[0][0] - 1])
                    new_interval.append(intervals[0])
                else:
                    new_interval.append(intervals[0])
            else:
                new_interval.append([intervals[intervals.index(i) - 1][1] + 1, i[0] - 1])
                new_interval.append(i)

            count += 1
        for i in samples:
            for j in new_interval:
                if i >= j[0] and i <= j[1]:
                    label.append(new_interval.index(j))

            if i > new_interval[len(new_interval) - 1][1]:
                label.append(len(new_interval))

        label_interval = []
        for elem in intervals:
            label_interval.append(new_interval.index(elem))

        process = list(set(label))
        result = []
        for g in label_interval:
            result.append(process.index(g))

        sample = []
        for i in samples:
            sample.append([i])

        samples = np.array(sample)

        if len(list(set(label))) > 1:
            labels = np.array(label)
            sample_silhouette_values = silhouette_samples(samples, labels)

            silhouette_per_cluster = {}
            for label in np.unique(labels):
                silhouette_per_cluster[label] = np.mean(sample_silhouette_values[labels == label])

            for i in label_interval:
                result_file.write(gene + '\t' + str(new_interval[i][0]) + ',' + str(new_interval[i][1]) + '\t' + str(
                    silhouette_per_cluster[i]) + '\n')
                ###If there is only one cluster in the transcript, the contour coefficient is 0.4 to 0.6
        elif len(list(set(label))) == 1 and len(new_interval) == 1:
            result_file.write(gene + '\t' + str(new_interval[0][0]) + ',' + str(new_interval[0][1]) + '\t' + str(
                format(random.uniform(0.4, 0.6), '.8f')) + '\n')
        else:
            result_file.write(gene + '\t' + str(new_interval[1][0]) + ',' + str(new_interval[1][1]) + '\t' + str(
                format(random.uniform(0.4, 0.6), '.8f')) + '\n')

    result_file.close()
    del df
    del all_variant

    ###Write the contour coefficient result into the interval result
    file_input = open(
        here + 'tmp/hot_cold_spot_modify(no-mutation-ratio)(2-sigama).txt', 'r')
    file_profile = open(
        here + 'tmp/profile-coefficient(no-mutation-ratio)(2-sigama).txt', 'r')
    file_result = open(
        here + 'tmp/hot_cold_result-and-profile-coefficient(no-mutation-ratio)(2-sigama).txt',
        'w')

    count = 0
    for line in file_input:
        flag = 0
        file_profile.seek(0)
        if count == 0:
            file_result.write(line.replace('\n', '') + '\t' + 'profile-coefficient' + '\n')
        for pro in file_profile:
            data = line.split('\t')
            coe = pro.split('\t')
            if data[0] != 'type':
                if data[2] == coe[0] and int(data[4]) == int(coe[1].split(',')[0]) and int(data[5]) == int(
                        coe[1].split(',')[1]):
                    file_result.write(line.replace('\n', '') + '\t' + coe[2])
                    flag = 1
                    break

        count += 1

    file_input.close()
    file_profile.close()
    file_result.close()

    print('ok')
    
if __name__ == '__main__':
    compute_profile_coefficient()
       
        


