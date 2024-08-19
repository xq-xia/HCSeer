#!/usr/bin/python3

### Made by xiaxingquan
### April 2024

# -*- coding: utf-8 -*-

import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from main import KDE

np.seterr(divide='ignore', invalid='ignore')

from scipy.stats import norm

np.random.seed(0)  # for reproducibility

here = '/data/xiaxq/topic_PM1/topic_PM1_code/database/'

def e_step(x, mean, std, w):
    """E-Step: Compute the responsibilities
    """
    # Compute the densities of the points under the two normal distributions
    prob = []
    for i in range(0, len(mean)):
        prob.append(norm(mean[i], std[i]).pdf(x) * w[i])

    # Normalize the probabilities
    prob_sum = 0.0
    for i in range(0, len(mean)):
        prob_sum += prob[i]
    for i in range(0, len(mean)):
        prob[i] = prob[i] / prob_sum

    return prob


def m_step(x, prob):
    """M-Step: Update the GMM parameters
    """
    # Update means
    mean = []
    for i in range(0, len(prob)):
        mean.append(np.dot(prob[i], x) / np.sum(prob[i]))

    # Update standard deviations
    std = []
    for i in range(0, len(prob)):
        std.append(np.sqrt(np.dot(prob[i], (x - mean[i]) ** 2) / np.sum(prob[i])))

    # Update mixing weights
    w = []
    for i in range(0, len(prob)):
        w.append(np.sum(prob[i]) / len(x))

    return mean, std, w


def gmm_em(x, max_iter=100, mean=[], std=[], w=[]):
    """Gaussian mixture model estimation using Expectation-Maximization
    """

    for i in range(max_iter):
        flag_mean = []
        flag_std = []
        flag_w = []
        # remove not exist mean,std and w
        for g in range(0, len(mean)):
            if std[g] != 0.0:
                flag_mean.append(mean[g])
                flag_std.append(std[g])
                flag_w.append(w[g])
        mean = flag_mean
        std = flag_std
        w = flag_w
        if len(mean) != 0 and 'nan' != str(mean[0]):
            prob = e_step(x, mean, std, w)
            mean, std, w = m_step(x, prob)
        else:
            return [], [], []

    return mean, std, w


##Interval overlap processing
def merge_intervals(intervals):
    # Sort intervals by their starting position
    intervals.sort(key=lambda x: x[0])

    # Initialize the merged interval list
    merged = []

    for interval in intervals:
        # If the merged list is empty or the current interval does not overlap with the last interval in the list
        if not merged or merged[-1][1] < interval[0]:
            merged.append(interval)
        else:
            # Otherwise, merge intervals
            if merged[-1][1] > interval[1]:
                merged[-1][1] = merged[-1][1]
            else:
                merged[-1][1] = interval[1]

    return merged


### result file

def EM():
    result_file = open(
        here + 'update_database/Initial_results_modify(1-sigama)_20240803', 'w')

    # get data
    df = pd.read_csv(here + 'update_database/final_variant.txt', sep='\t',
                     low_memory=False)
    genename = df.iloc[:, 6]
    print(len(genename.unique()))
    P = df[df['type'] == 'PLP']
    B = df[df['type'] == 'BLB']

    print('PLP', P['genename'].unique())
    print('BLB', B['genename'].unique())

    for elem in genename.unique():
        gene = elem
        flag = df[df['genename'] == gene]
        chromosome = flag['chr'].unique()[0]
        iso = flag['iso'].unique()[0]
        max_position = np.array(flag['aa_position']).max()

        ##variant
        PLP = flag[flag['type'] == 'PLP']
        BLB = flag[flag['type'] == 'BLB']

        ##variant position
        PLP_pos = np.array(PLP['aa_position'])
        BLB_pos = np.array(BLB['aa_position'])

        #######analysis PLP variants

        PLP_primary_data = []
        # # print(len(flag))
        for i in PLP_pos:
            PLP_primary_data.append(i)
        # print('primary',PLP_primary_data)
        PLP_x_all = np.linspace(0, max_position,
                                max_position)

        ### prevent variant length from being 0, resulting in a denominator of 0
        if len(PLP_primary_data) > 0:
            kde = KDE(PLP_primary_data)
            #    print('Optimal bandwidth:', kde.hopt)
            PLP_density_estimation_all = np.array([kde.compute_density(xi) for xi in PLP_x_all])
        else:
            PLP_density_estimation_all = np.array([0 for xi in PLP_x_all])

        # ##Process the original data to remove locations that clearly do not constitute hotspots

        data_splice = []
        flag = 0
        for i in range(len(PLP_primary_data) - 1):
            if PLP_primary_data[i + 1] - PLP_primary_data[i] >= 200:
                data_splice.append(PLP_primary_data[flag:i + 1])
                flag = i + 1
        data_splice.append(PLP_primary_data[flag:])
        data = []
        for i in data_splice:
            if len(i) > 1:
                data += i
        x = np.linspace(0, max_position,
                        max_position)
        if len(data) != 0:
            kde = KDE(data)
            #   print('Optimal bandwidth:', kde.hopt)
            density_estimation = np.array([kde.compute_density(xi) for xi in x])
            ###Obtain maximum kernel density

            max = []
            if density_estimation[1] < density_estimation[0]:
                max.append(0)
            for i in range(1, len(density_estimation) - 1):
                if density_estimation[i] > density_estimation[i - 1] and density_estimation[i] > density_estimation[
                    i + 1]:
                    max.append(i)
            if density_estimation[len(density_estimation) - 2] < density_estimation[len(density_estimation) - 1]:
                max.append(len(density_estimation) - 1)

            ###Obtain initial mean, standard deviation, and weight

            mean = max
            std = []
            w = []
            for elem in max:
                left, right = 0, max[len(max) - 1]
                if elem != 0:
                    for i in range(elem, -1, -1):
                        if i == 0:
                            left = 0
                            break
                        if density_estimation[i - 1] >= density_estimation[i]:
                            left = i
                            break
                if elem != max[len(max) - 1]:
                    for i in range(elem, max[len(max) - 1] + 1, 1):
                        if i == max[len(max) - 1]:
                            right = max[len(max) - 1]
                            break
                        if density_estimation[i + 1] >= density_estimation[i]:
                            right = i
                            break
                # print('left', left, 'right', right)
                std.append((right - left) / 2)
                w.append(1 / len(max))
            #
            ###Perform expected value maximization iteration

            final_dist_params = gmm_em(data, max_iter=30, mean=mean, std=std, w=w)
            PLP_result = []
            for i in range(len(final_dist_params[0])):
                if math.isnan(final_dist_params[0][i]):
                    final_dist_params[0][i] = 0
                if math.isnan(final_dist_params[1][i]):
                    final_dist_params[1][i] = 0
                flag = [round(final_dist_params[0][i] - math.ceil(final_dist_params[1][i])) - 1,
                        round(final_dist_params[0][i] + math.ceil(final_dist_params[1][i])) + 1]
                if flag[0] < 0:
                    flag[0] = 0
                if flag[1] > max_position:
                    flag[1] = max_position
                PLP_result.append(flag)
            PLP_result = merge_intervals(PLP_result)

        #    print('hotspot', PLP_result)
        #    print(data)
        else:
            density_estimation = np.array([0 for xi in x])
            PLP_result = []

        ### get hotspot area
        PLP_d = []
        for i in range(len(PLP_result)):
            j = PLP_result[i][0]
            while j <= PLP_result[i][1]:
                PLP_d.append(j)
                j += 0.1

        ### assignment
        PLP_y = []
        for i in range(len(PLP_d)):
            PLP_y.append(-0.0025)

        ####### similarly for BLB

        BLB_primary_data = []
        # # print(len(flag))
        for i in BLB_pos:
            BLB_primary_data.append(i)
        # print('primary', BLB_primary_data)
        BLB_x_all = np.linspace(0, max_position,
                                max_position)

        ### prevent variant length from being 0, resulting in a denominator of 0
        if len(BLB_primary_data) > 0:
            kde = KDE(BLB_primary_data)
            #    print('Optimal bandwidth:', kde.hopt)
            BLB_density_estimation_all = np.array([kde.compute_density(xi) for xi in BLB_x_all])
        else:
            BLB_density_estimation_all = np.array([0 for xi in BLB_x_all])

        # ##Process the original data to remove locations that clearly do not constitute hotspots

        data_splice = []
        flag = 0
        for i in range(len(BLB_primary_data) - 1):
            if BLB_primary_data[i + 1] - BLB_primary_data[i] >= 200:
                data_splice.append(BLB_primary_data[flag:i + 1])
                flag = i + 1
        data_splice.append(BLB_primary_data[flag:])
        data = []
        for i in data_splice:
            if len(i) > 1:
                data += i
        # print('splice', data)
        x = np.linspace(0, max_position,
                        max_position)
        if len(data) != 0:
            kde = KDE(data)
            #    print('Optimal bandwidth:', kde.hopt)
            density_estimation = np.array([kde.compute_density(xi) for xi in x])
            ###Obtain maximum kernel density

            max = []
            if density_estimation[1] < density_estimation[0]:
                max.append(0)
            for i in range(1, len(density_estimation) - 1):
                if density_estimation[i] > density_estimation[i - 1] and density_estimation[i] > density_estimation[
                    i + 1]:
                    max.append(i)
            if density_estimation[len(density_estimation) - 2] < density_estimation[len(density_estimation) - 1]:
                max.append(len(density_estimation) - 1)

            ###Obtain initial mean, standard deviation, and weight

            mean = max
            std = []
            w = []
            for elem in max:
                left, right = 0, max[len(max) - 1]
                if elem != 0:
                    for i in range(elem, -1, -1):
                        if i == 0:
                            left = 0
                            break
                        if density_estimation[i - 1] >= density_estimation[i]:
                            left = i
                            break
                if elem != max[len(max) - 1]:
                    for i in range(elem, max[len(max) - 1] + 1, 1):
                        if i == max[len(max) - 1]:
                            right = max[len(max) - 1]
                            break
                        if density_estimation[i + 1] >= density_estimation[i]:
                            right = i
                            break
                #       print('left', left, 'right', right)
                std.append((right - left) / 2)
                w.append(1 / len(max))
            #
            ###Perform expected value maximization iteration

            final_dist_params = gmm_em(data, max_iter=30, mean=mean, std=std, w=w)
            BLB_result = []
            for i in range(len(final_dist_params[0])):
                if math.isnan(final_dist_params[0][i]):
                    final_dist_params[0][i] = 0
                if math.isnan(final_dist_params[1][i]):
                    final_dist_params[1][i] = 0
                flag = [round(final_dist_params[0][i] - math.ceil(final_dist_params[1][i])) - 1,
                        round(final_dist_params[0][i] + math.ceil(final_dist_params[1][i])) + 1]
                if flag[0] < 0:
                    flag[0] = 0
                if flag[1] > max_position:
                    flag[1] = max_position
                BLB_result.append(flag)
            BLB_result = merge_intervals(BLB_result)

        #   print('coldspot', BLB_result)
        #   print(data)
        else:
            density_estimation = np.array([0 for xi in x])
            BLB_result = []

        ### get coldspot area
        BLB_d = []
        for i in range(len(BLB_result)):
            j = BLB_result[i][0]
            while j <= BLB_result[i][1]:
                BLB_d.append(j)
                j += 0.1

        ### assignment
        BLB_y = []
        for i in range(len(BLB_d)):
            BLB_y.append(-0.001)

        #### cumpute variant number in every area
        ###hotspot
        for elem in PLP_result:
            p_count = 0
            b_count = 0
            for plp in PLP_primary_data:
                if plp >= elem[0] and plp <= elem[1]:
                    p_count += 1

            for blb in BLB_primary_data:
                if blb >= elem[0] and blb <= elem[1]:
                    b_count += 1
            elem.append(p_count)
            elem.append(b_count)
            ratio = 1 - (b_count + 1 / 2) / (p_count + 1 / 2)
            elem.append(round(ratio, 3))

        ###coldspot
        for elem in BLB_result:
            p_count = 0
            b_count = 0
            for plp in PLP_primary_data:
                if plp >= elem[0] and plp <= elem[1]:
                    p_count += 1

            for blb in BLB_primary_data:
                if blb >= elem[0] and blb <= elem[1]:
                    b_count += 1
            elem.append(p_count)
            elem.append(b_count)
            ratio = 1 - (p_count + 1 / 2) / (b_count + 1 / 2)
            elem.append(round(ratio, 3))

        PLP_txt = '(PLP)'
        for elem in PLP_result:
            txt = ':'.join(str(word) for word in elem)
            PLP_txt += txt + ';'

        BLB_txt = '(BLB)'
        for elem in BLB_result:
            txt = ':'.join(str(word) for word in elem)
            BLB_txt += txt + ';'

        result_file.write(chromosome + '\t' + gene + '\t' + iso + '\t' + PLP_txt + '\t' + BLB_txt + '\n')

    result_file.close()

    '''
    ### visualization of results
    line1, = plt.plot(PLP_x_all, PLP_density_estimation_all, color='r', marker='.', linestyle='-', markersize=2, alpha=0.5, linewidth=3)
    line2, = plt.plot(PLP_d, PLP_y, color='black', marker='s', linestyle='-', markersize=8, alpha=0.5, linewidth=8)

    line3, = plt.plot(BLB_x_all, BLB_density_estimation_all, color='blue', marker='.', linestyle='-', markersize=2, alpha=0.5, linewidth=3)
    line4, = plt.plot(BLB_d, BLB_y, color='green', marker='s', linestyle='-', markersize=8, alpha=0.5, linewidth=8)

    plt.xlim(0,max_position)
    #hopspot area
    plt.legend(['PLP Nuclear density','hotspot area','BLB Nuclear density','coldspot area'])
    title = 'Nuclear density for ' + gene
    plt.title(title,fontsize=18, color='black')
    plt.xlabel(u'Amino acid position', fontsize=14, color='r')
    plt.ylabel(u'Nuclear density', fontsize=14, color='b')
    plt.show()
    '''

if __name__ == '__main__':
    EM()