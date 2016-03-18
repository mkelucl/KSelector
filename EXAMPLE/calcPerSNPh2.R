# Must set this
print("Started R script to generate per SNP h2 contributions")

# _________________________________________________________

# set to current working dir
currentDir = '/home/XXX';
setwd(currentDir)


## _________________________________
# Load in data, total LDAK h2, Weights and per SNP effects
## _________________________________


# load iteration:
iteration=""
if(file.exists('num')) {
iteration = as.numeric( read.table('num',header=F) )
}

# load in total LDAK h2 from my java processed file
h2File = read.table(paste('h2_',iteration,'.txt',sep=""),header=F)

sectionsname ='sections'

# get total h2 we are looking for
h2 = as.character(h2File$V1[length(h2File$V1)])
myList=unlist(strsplit(h2, split='"'))
LDAK_totalh2 = as.numeric(myList[2]) # 
print( paste("LDAK non BG Region h2 is: " , LDAK_totalh2, sep=""))


# get the weights
SNP_Weights = read.table(paste(sectionsname,iteration,'/weightsALL',sep=""),header=F);  # this holds the SNP weights

# get per SNP effects
ldak_snpreg_effects_table = read.table(paste('snpfolder',iteration,'/regress1',sep=""),header=T);  # this is the no weight set
ldak_snpreg_effects = ldak_snpreg_effects_table$REML_Her;  # this holds the SNP h2 contributions for each SNP



# _________________________________
# estimate the top SNP h2s from the per SNP LDAK effects
# _________________________________

# Weight the per SNP effects by LD
ldak_snpreg_effects_weighted = ldak_snpreg_effects * SNP_Weights$V1
diff = sum(ldak_snpreg_effects_weighted) - LDAK_totalh2
ldak_snpreg_scaled_weighted2 = ldak_snpreg_effects_weighted - (diff / length(ldak_snpreg_effects_weighted))

snph2ContributionFrame= data.frame(SNP_Weights$V5,ldak_snpreg_scaled_weighted2)

filename="SNPH2s.csv"
write.table(snph2ContributionFrame, filename, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
print( paste("written " , length(ldak_snpreg_scaled_weighted2), " predictors to file: ", filename, sep=""))


