#!/bin/bash

### Made by xiaxingquan
### May 2024

###The purpose of this file is to convert the variation of hg19 version to hg38 version when commenting on the variation of hg19 version

here=/data/xiaxq/topic_PM1/topic_PM1_code

## params 0: script name 1:input file 2:output file 3:convert type
para_0=$0
para_1=$1
para_2=$2
para_3=$3

### convert bed file format
awk -F '\t'  '{if(NF == 5 ){print $1FS$2-1FS$3FS$4FS$5}else{print $1FS$2-1FS$3FS$4FS$5FS$6}}' $para_1 > $here/database/tmp/hg19-to-hg38-middle-file.bed

### run crossMap script
CrossMap  bed   $here/database/chain_files/$3.gz  $here/database/tmp/hg19-to-hg38-middle-file.bed  $here/database/tmp/hg19-to-hg38-middle-file2.bed

### convert vcf file format
awk -F '\t' '{if(NF == 5 ){print $1FS$2+1FS$3FS$4FS$5}else{print $1FS$2+1FS$3FS$4FS$5FS$6}}'  $here/database/tmp/hg19-to-hg38-middle-file2.bed  > $para_2

### Remove the generated intermediate files
rm  $here/database/tmp/hg19-to-hg38-middle-file.bed
rm  $here/database/tmp/hg19-to-hg38-middle-file2.bed
rm  $here/database/tmp/hg19-to-hg38-middle-file2.bed.unmap