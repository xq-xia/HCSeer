# 变异热点冷点区域

#### 介绍
计算变异热点区域及冷点区域


####需要安装的python包
1、numpy (>=v1.26.4)
2、pandas (>=v2.2.2)
3、scikit-learn (v>=1.4.2)
4、scipy (v>=1.13.1)

#### 安装教程

这些脚本不需要编译。

#### 使用说明

- 00_database-processing.sh
  对clinvar数据库的数据进行预处理
- 01_clinvar-annotation.sh
  对预处理后的数据用annovar进行注释
- 03_annotation-data-processing.py
  对注释生成的数据进行处理
- 04_main.py
  核密度估计算法
- 05_EM.py
  期望值最大算法
- 06_Initial-result-processing.py
  处理初始结果得到最终结果
- 07_compute_p-value.py
  计算p-value
- 08_profile-coefficient.py
  计算轮廓系数
- Auto-annotation.py
  注释脚本
- hg19-to-hg38.sh
  hg19基因组版本坐标转换为hg38版本的脚本
- inAutoPVS1.py
  与AutoPVS1比较
- PP2_BP1_PM1_PS1.py
  结合PP2、PM5、PS1、BP1
- process_function.py
  一些文件处理函数
- result-analysis.py
  结果分析
