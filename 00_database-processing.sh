#!/bin/bash

### Made by xiaxingquan
### April 2024

here=/data/xiaxq/topic_PM1/topic_PM1_code
batch="20240407"

#### ClinVar

# VCF with variants from ClinVar can be downloaded here: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
# It has be put in the following directory:
# mkdir -p $here/database

grep -P "CLNSIG=Pathogenic;|CLNSIG=Likely_pathogenic;|CLNSIG=Pathogenic/Likely_pathogenic;" $here/database/clinvar_$batch.vcf |grep -E "CLNREVSTAT=criteria_provided,_multiple_submitters|reviewed_by_expert_panel|CLNREVSTAT=practice_guideline" | awk -F"\t" '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 }' | awk -F"\t" '{split($6,a,";");for(j in a){ split(a[j],m,"=");if(m[1] == "GENEINFO"){split(m[2],g,":");print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" g[1];}}}' > $here/database/update_database/clinvar-$batch-PLP-temp.vcf

grep -P "CLNSIG=Likely_benign;|CLNSIG=Benign;|CLNSIG=Benign/Likely_benign;" $here/database/clinvar_$batch.vcf | grep -E "CLNREVSTAT=criteria_provided,_multiple_submitters|reviewed_by_expert_panel|CLNREVSTAT=practice_guideline" | awk -F"\t" '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 }' | awk -F"\t" '{split($6,a,";");for(j in a){ split(a[j],m,"=");if(m[1] == "GENEINFO"){split(m[2],g,":");print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" g[1];}}}' > $here/database/update_database/clinvar-$batch-BLB-temp.vcf

grep -P "CLNSIG=Uncertain_significance;" $here/database/clinvar_$batch.vcf | grep -E "CLNREVSTAT=criteria_provided,_multiple_submitters|reviewed_by_expert_panel|CLNREVSTAT=practice_guideline" | awk -F"\t" '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 }' | awk -F"\t" '{split($6,a,";");for(j in a){ split(a[j],m,"=");if(m[1] == "GENEINFO"){split(m[2],g,":");print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" g[1];}}}' > $here/database/update_database/clinvar-$batch-VUS-temp.vcf

grep -P "CLNSIG=Conflicting_classifications_of_pathogenicity;" $here/database/clinvar_$batch.vcf | grep -E "CLNREVSTAT=criteria_provided,_multiple_submitters|reviewed_by_expert_panel|CLNREVSTAT=practice_guideline" | awk -F"\t" '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 }' | awk -F"\t" '{split($6,a,";");for(j in a){ split(a[j],m,"=");if(m[1] == "GENEINFO"){split(m[2],g,":");print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" g[1];}}}' > $here/database/update_database/clinvar-$batch-CON-temp.vcf

cat  $here/database/head.txt  $here/database/update_database/clinvar-$batch-VUS-temp.vcf >$here/database/update_database/clinvar-$batch-VUS.vcf
cat  $here/database/head.txt  $here/database/update_database/clinvar-$batch-PLP-temp.vcf >$here/database/update_database/clinvar-$batch-PLP.vcf
cat  $here/database/head.txt  $here/database/update_database/clinvar-$batch-BLB-temp.vcf >$here/database/update_database/clinvar-$batch-BLB.vcf
cat  $here/database/head.txt  $here/database/update_database/clinvar-$batch-CON-temp.vcf >$here/database/update_database/clinvar-$batch-CON.vcf
rm $here/database/update_database/clinvar-$batch-PLP-temp.vcf
rm $here/database/update_database/clinvar-$batch-BLB-temp.vcf
rm $here/database/update_database/clinvar-$batch-VUS-temp.vcf
rm $here/database/update_database/clinvar-$batch-CON-temp.vcf

####MANE

#isoform from MANE can be downloaded here:https: https://ftp.ncbi.nlm.nih.gov/refseq/MANE/




####variant_summary

#variant_summary.txt can be downloaded here:https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
grep GRCh38 $here/database/variant_summary.txt | grep -i pathogenic | grep -v 'Conflicting classifications of pathogenicity' | cut -f3,5,26 > $here/database/variant_summary_PLP.txt
grep -i benign $here/database/variant_summary,txt |grep GRCh38 | cut -f3,5,26 > $here/database/variant_summary_BLB.txt
head -n 1 $here/database/variant_summary.txt | cut -f3,5,26 > $here/database/variant_summary_head.txt
cat $here/database/variant_summary_head.txt $here/database/variant_summary_PLP.txt $here/database/variant_summary_BLB.txt > $here/database/variant_summary_PLP_BLB.txt

rm $here/database/variant_summary_BLB.txt
rm $here/database/variant_summary_PLP.txt
rm $here/database/variant_summary_head.txt