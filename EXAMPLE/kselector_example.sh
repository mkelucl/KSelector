#!/bin/bash
#$ -S /bin/bash
#$ -o /home/XXX/test_ex.out
#$ -e /home/XXX/test_ex.err
#$ -l tmem=20G,h_vmem=20G
#$ -cwd
#$ -l h_rt=120:0:0
#$ -V

##############################################################################################################################
##############################################################################################################################
# EDIT THIS TO YOUR LOCATION: 
# ( YOU WILL ALSO NEED TO EDIT THE  LOCATIONS IN THE FILE:'calcPerSNPh2.R' )

PATH=~/bin:$PATH
PATH=/share/apps/R/bin:$PATH
java=/share/apps/jdk1.8.0_25/jre/bin/java

homeBase='/home/XXX/'
kselector=$homeBase$'KSelector.jar'
ldak2016=$homeBase$'ldak2016.5.94' 
plink=$homeBase$'plink'


### DONT FORGET TO  FILTER FOR DUPLICATE SNPS, OTHERWISE THE WEIGHTS WILL BE PERMUTED!  ###
# ./plink --bfile ORIG_DATASET --bp-space 1 --make-bed --out ORIG_QC

# THESE SHOULD POINT TO THE TRAINING/TEST SETS AND THE PHENOTYPES
trainingSet='YOUR_TRAINING_PLINK'
testSet='YOUR_TEST_PLINK'
phenotypes='YOURPHENO.pheno'

##############################################################################################################################
##############################################################################################################################


##################### 
echo 'STAGE 1: AMBLUP gross filter '
#####################
date

arguments=' --cut-genes chunks_amb --chunks-bp 75000 --bfile '$trainingSet$' --ignore-weights YES'
$ldak2016 $arguments

arguments=' --calc-genes-reml chunks_amb --bfile '$trainingSet$' --pheno '$phenotypes$' --ignore-weights YES --partition 1'
$ldak2016 $arguments

arguments=' --join-genes-reml chunks_amb --bfile '$trainingSet$' --sig1 1e-5 --sig2 0.01'
$ldak2016 $arguments

arguments=' --calc-kins-direct chunks_amb/region0 --bfile '$trainingSet$' --extract chunks_amb/region0 --ignore-weights YES'
$ldak2016 $arguments

file="chunks_amb/region_number.txt"
read -d $'\x04' regionnumber < "$file"

arguments=' --reml amblup --grm chunks_amb/region0 --pheno '$phenotypes$' --region-number '$regionnumber$' --region-prefix chunks_amb/region --bfile '$trainingSet$' --ignore-weights YES'
$ldak2016 $arguments

arguments=' --calc-blups amblup --grm chunks_amb/region0 --remlfile amblup.reml --bfile '$trainingSet$''
$ldak2016 $arguments

# prediction for basic AMBLUP
arguments=' --calc-scores amblup --scorefile amblup.blup --bfile '$testSet$''
$ldak2016 $arguments


# merge all AMBLUP surviving regions into single SNP extract to be used to obtain the weights/per SNP h2s
# parse the LDAK h2 into an easier accessible text file, this is the H2 that the regions are supposed to have it total when calculated individually by amblup (this will be GREATER than if we just run it on the full list)
arguments=' *parseh2 amblup.reml h2.txt'
$java -Xms2G -Xmx2G -jar $kselector $arguments
mv h2.txt 'h2_.txt'

# save backup copy of regions, as KSelector will overwrite the original regions
cp -rf chunks_amb 'chunks_amb_bak'

# this produces the SNP extract list -> MBLUPSNPs.csv
arguments=' *getregsnps chunks_amb/region'
$java -Xms2G -Xmx2G -jar $kselector $arguments

arguments=' --bfile '$trainingSet$' --extract MBLUPSNPs.csv --make-bed --out train_amblupregions'
$plink $arguments


##################### 
echo 'STAGE 2: calculate WEIGHTS and per SNP h2s - only on the AMBLUP surviving regions'
#####################

#### calculate weights
arguments=' --cut-weights sections --bfile train_amblupregions'
$ldak2016 $arguments

file='sections/section_number'
read -d $'\x04' numsections < "$file"

for k in `seq 1 $((numsections))`; do
# go through all sections
echo 'calculating weights for section'$k

arguments=' --calc-weights sections --bfile train_amblupregions --section '$k
$ldak2016 $arguments

done 

arguments=' --join-weights sections'
$ldak2016 $arguments


#### calculate per SNP h2s (unadjusted)
arguments=' --cut-genes snpfolder --chunks-bp 1 --bfile train_amblupregions --ignore-weights YES'
$ldak2016 $arguments

arguments=' --calc-genes-reml snpfolder --bfile train_amblupregions --pheno '$phenotypes$' --grm chunks_amb/region0 --ignore-weights YES --partition 1'
$ldak2016 $arguments


# get per SNP h2 contributions
# load in current h2 file 
source 'h2_.txt'

echo 'ALL H2 from non-bg regions is: '$H2_REGIONS

# R: this relies on sections/weightsALL, snpfolder/regress1 being ready -> produces a SNPH2s.csv, that has each SNPs estimated h2 contribution (original order)
Rscript calcPerSNPh2.R


##################### 
echo 'STAGE 3: Fine filter via KSelector: h2 guided regional lasso '
#####################

# shrinkageFactor='0.7'
# arguments=' *inputplink train_amblupregions *inputphenotype phenotypes5_train.pheno *snp2h2 SNPH2s.csv *h2helper amblup.reml chunks_amb/region *h2 '$H2_REGIONS$' *shfact '$shrinkageFactor
arguments=' *inputplink train_amblupregions *inputphenotype phenotypes5_train.pheno *snp2h2 SNPH2s.csv *h2helper amblup.reml chunks_amb/region *h2 '$H2_REGIONS
$java -Xms12G -Xmx12G -jar $kselector $arguments
# produces -> lassoresults.csv, lassoresults.csv.log, lassoresults.csv.extract


##################### 
echo '6. LDAK: perform AMBLUP on the KSelect thinned regions'
#####################

file="chunks_amb/region_number.txt"
read -d $'\x04' regionnumber < "$file"

arguments=' --reml amblup_ks --grm chunks_amb/region0 --pheno '$phenotypes$' --region-number '$regionnumber$' --region-prefix chunks_amb/region --bfile '$trainingSet$' --ignore-weights YES'
$ldak2016 $arguments

arguments=' --calc-blups amblup_ks --grm chunks_amb/region0 --remlfile amblup_ks.reml --bfile '$trainingSet$''
$ldak2016 $arguments

# prediction for KS-BLUP
arguments=' --calc-scores amblup_ks --scorefile amblup_ks.blup --bfile '$testSet$''
$ldak2016 $arguments



## your profile scores will be in these 2 files:
# amblup_ks.profile
# amblup.profile
date
