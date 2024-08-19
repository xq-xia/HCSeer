import pandas as pd
import numpy as np

import openpyxl # 导入openpyxl模块

def penetrance_extract():
    wb = openpyxl.load_workbook('C:/Users/xiaxq/Desktop/topic-PM1-tep-file/2-final20240815.xlsx') # 创建workbook对象
    ws = wb.active# 得到当前活跃表单的对象（打开该xlsx文件，直接看到的表单就为活跃表单）

    result = open('C:/Users/xiaxq/Desktop/topic-PM1-tep-file/penetrance.txt','w')

    g = []
    print('ok')

    count = 0
    for line in ws:
        result.write(str(line[0].value) + '\t' + str(line[1].value) + '\t' + str(line[2].value) + '\t' + str(line[3].value) + '\t' + str(line[4].value) + '\t' + str(line[12].value) + '|' + str(line[13].value) + '|' + str(line[14].value) + '|' + str(line[15].value) + '|' + str(line[16].value) + '\n' )

    print('ok')
    result.close()

'''
chr	start	end	alt	ref	otherInfo	hot_cold_score
4	1799344	1799344	G	A	0|4|17|43374|0	0
'''
def result_process():
    input = open('C:/Users/xiaxq/Desktop/topic-PM1-tep-file/penetrance_annotation.txt', 'r')
    result = open('C:/Users/xiaxq/Desktop/topic-PM1-tep-file/penetrance_annotation_presscoss.txt', 'w')
    count = 0
    for line in input:
        if count != 0:
            data = line.split('\t')[-2].split('|')
            if data[0] == 'NA' or data[1] == 'NA' or data[2] == 'NA' or data[3] == 'NA':
                continue
            if float(data[0]) < 0 or float(data[1]) < 0 or float(data[2]) < 0 or float(data[3]) < 0:
                continue
            if float(data[0]) > float(data[1]) or float(data[2]) > float(data[3]):
                continue
            result.write(line)
        count += 1
    input.close()
    result.close()


if __name__ == '__main__':
    result_process()