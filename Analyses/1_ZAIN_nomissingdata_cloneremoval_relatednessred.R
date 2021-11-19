#############This is the first script for the Zamia integrifolia in situ/ex situ project
############We are comparing ex situ collection of wild diversity 
###########For the rare cycad Zamia integrifolia
##########This script will focus on individual reduction based on missing data and relatedness

###################################
########## Libraries ##############
###################################
##load in packages
library(adegenet)
library(diveRsity)
library(poppr)
library(hierfstat)

##################################
####### Load in data files #######
##################################
##setting your working directory 
ZAIN_wd <- "G:\\Shared drives\\Emily_Schumacher\\ZAIN"

##set working directory
setwd(ZAIN_wd)

##convert to a genind file 
arp2gen("\\Data_Files\\ZAIN_garden_wild\\ZAIN_garden_wild.arp")

##load in genind 
ZAIN_garden_wild_gen <- read.genepop("Data_Files\\ZAIN_garden_wild\\ZAIN_garden_wild.gen", ncode = 3)

##ZAIN score df 
ZAIN_garden_wild_df <- read.csv(paste0(ZAIN_wd, "\\Data_Files\\ZAIN_garden_wild\\ZAIN_garden_wild_df.csv"))

##organize genind file 
levels(ZAIN_garden_wild_gen@pop) <- unique(ZAIN_garden_wild_df$Pop)

##name the individuals in the genind 
rownames(ZAIN_garden_wild_gen@tab) <- ZAIN_garden_wild_df$Sample.Name

#########################################################
##### Remove clones, missing data, and siblings #########
#########################################################
##remove individuals with too much missing data 
ZAIN_garden_wild_gen <- missingno(ZAIN_garden_wild_gen, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE)

#First put into poppr format
ZAIN_garden_clones <- as.genclone(ZAIN_garden_gen)
strata(ZAIN_garden_clones) <- other(ZAIN_garden_clones)$population_hierarchy[-1]
list_a<-mlg.id(ZAIN_garden_clones)
#Function to pull out individual indices where clone length greater than 1
clone_index<-which(sapply(list_a,function(x) length(x)>1))
write.csv(list_a[clone_index],file=paste0(species_names[sp],"clones.csv"))

#This removes clones and then saves as new file for Genealex if desired
popr_nocl<-clonecorrect(ZAIN_garden_clones,strata=~Pop)
#genind2genalex(genclone2genind(popr_nocl),file="QH_clone_free.csv")

#Create genpop and genind objects that now have no clones- GI_nocl, GP_nocl
ZAIN_garden_noclones <-genclone2genind(popr_nocl)

##create poppr file 
ZAIN_poppr <- poppr(ZAIN_garden_noclones)

##expected heterozygosity 
ZAIN_hexp <- ZAIN_poppr[1:10, 10]

##
ZAIN_all_rich <- colMeans(allelic.richness(ZAIN_garden_noclones)$Ar)
