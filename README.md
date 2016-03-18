KSelector
===============

A java based application that applies a fine filter to the AMBLUP surviving regions in order to reduce the model size to the minimum number of SNPs that still maintain the same predictive power as the AMBLUP superset.


This is a command line utility that requires the following parameters:

1. the input PLINK binary datasets:
*inputplink YOURTRAINGENOTYPES

2. the input phenotype data :
*inputphenotype YOURTRAINGENOTYPES

3. the per SNP heritability contributions:
*snp2h2 SNPH2s.csv

4. the location of the regions that AMBLUP found:
*h2helper chunks_amb/region

5. the overall heritability in the dataset:
*h2 0.5


Full command line:
java -jar KSelector.jar *inputplink YOURTRAINGENOTYPES *inputphenotype YOURTRAININGPHENOS.pheno *snp2h2 SNPH2s.csv *h2helper amblup.reml chunks_amb/region *h2 0.5


PS:
an example of how the application could be run from a cluster can be found under:

EXAMPLE/kselector_example.sh
