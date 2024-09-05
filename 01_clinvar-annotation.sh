#!/bin/bash

### Made by Xiaxingquan
### April 2024

here=/data/xiaxq
batch=20240407

###############################################
##Download the Human Gene Reference Database

perl $here/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene $here/annovar/humandb/

###database processing
sed -i '1d' $here/topic_PM1/topic_PM1_code/database/clinvar-$batch-BLB.vcf
cut -f1,2,4,5 $here/topic_PM1/topic_PM1_code/database/clinvar-$batch-BLB.vcf | awk -F "\t" '{m=0;len=length($3);m=$2+len-1;chr=gsub(/chr/, $1,"");print $1 "\t" $2 "\t" m "\t" $3 "\t" $4}' > $here/topic_PM1/topic_PM1_code/database/BLB.vcf
sed -i '1d' $here/topic_PM1/topic_PM1_code/database/clinvar-$batch-PLP.vcf
cut -f1,2,4,5 $here/topic_PM1/topic_PM1_code/database/clinvar-$batch-PLP.vcf | awk -F "\t" '{m=0;len=length($3);m=$2+len-1;chr=gsub(/chr/, $1,"");print $1 "\t" $2 "\t" m "\t" $3 "\t" $4}' > $here/topic_PM1/topic_PM1_code/database/PLP.vcf
sed -i '1d' $here/topic_PM1/topic_PM1_code/database/clinvar-$batch-VUS.vcf
cut -f1,2,4,5 $here/topic_PM1/topic_PM1_code/database/clinvar-$batch-VUS.vcf | awk -F "\t" '{m=0;len=length($3);m=$2+len-1;chr=gsub(/chr/, $1,"");print $1 "\t" $2 "\t" m "\t" $3 "\t" $4}' > $here/topic_PM1/topic_PM1_code/database/VUS.vcf
sed -i '1d' $here/topic_PM1/topic_PM1_code/database/clinvar-$batch-CON.vcf
cut -f1,2,4,5 $here/topic_PM1/topic_PM1_code/database/clinvar-$batch-CON.vcf | awk -F "\t" '{m=0;len=length($3);m=$2+len-1;chr=gsub(/chr/, $1,"");print $1 "\t" $2 "\t" m "\t" $3 "\t" $4}' > $here/topic_PM1/topic_PM1_code/database/CON.vcf

### annotation by annovar
perl 	/data/xiaxq/annovar/annotate_variation.pl -geneanno -buildver hg38 -dbtype refGene -outfile /data/xiaxq/topic_PM1/topic_PM1_code/database/clinvar_annotation/annotation-PLP.refGene -exonsort -nofirstcodondel /data/xiaxq/topic_PM1/topic_PM1_code/database/PLP.vcf /data/xiaxq/annovar/humandb/

perl 	/data/xiaxq/annovar/annotate_variation.pl -geneanno -buildver hg38 -dbtype refGene -outfile /data/xiaxq/topic_PM1/topic_PM1_code/database/clinvar_annotation/annotation-BLB.refGene -exonsort -nofirstcodondel /data/xiaxq/topic_PM1/topic_PM1_code/database/BLB.vcf /data/xiaxq/annovar/humandb/

perl 	/data/xiaxq/annovar/annotate_variation.pl -geneanno -buildver hg38 -dbtype refGene -outfile /data/xiaxq/topic_PM1/topic_PM1_code/database/clinvar_annotation/annotation-VUS.refGene -exonsort -nofirstcodondel /data/xiaxq/topic_PM1/topic_PM1_code/database/VUS.vcf /data/xiaxq/annovar/humandb/

perl 	/data/xiaxq/annovar/annotate_variation.pl -geneanno -buildver hg38 -dbtype refGene -outfile /data/xiaxq/topic_PM1/topic_PM1_code/database/clinvar_annotation/annotation-CON.refGene -exonsort -nofirstcodondel /data/xiaxq/topic_PM1/topic_PM1_code/database/CON.vcf /data/xiaxq/annovar/humandb/


###remove unnecessary files
rm $here/topic_PM1/topic_PM1_code/database/PLP.vcf
rm $here/topic_PM1/topic_PM1_code/database/BLB.vcf
rm $here/topic_PM1/topic_PM1_code/database/VUS.vcf
rm $here/topic_PM1/topic_PM1_code/database/CON.vcf
