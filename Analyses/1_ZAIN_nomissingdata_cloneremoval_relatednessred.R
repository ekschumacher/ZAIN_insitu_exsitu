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
library(Demerelate)

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

#############################################################
######## Remove individuals with too much missing data ######
#############################################################
##remove individuals with too much missing data 
ZAIN_garden_wild_gen <- missingno(ZAIN_garden_wild_gen, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE)

##################################
########## Remove clones #########
##################################
##Create a list of clones from the ZAIN inds 
ZAIN_clones <- as.genclone(ZAIN_garden_wild_gen)
##setting the population levels for the genind   
strata(ZAIN_clones) <- other(ZAIN_clones)$population_hierarchy[-1]
##create MLG list 
ZAIN_mlg_list <- mlg.id(ZAIN_clones)
#Function to pull out individual indices where clone length greater than 1
ZAIN_clone_index <- which(sapply(ZAIN_mlg_list, function(x) length(x)>1))
write.csv(ZAIN_mlg_list[ZAIN_clone_index],file = "C:\\Users\\eschumacher\\Documents\\GitHub\\ZAIN_insitu_exsitu\\Data_Files\\ZAIN_data_frames\\ZAIN_clone_list.csv")

#This removes clones and then saves as new file for Genlex if desired
ZAIN_nocl <- clonecorrect(ZAIN_garden_wild_gen, strata=~Pop)
genind2genalex(genclone2genind(ZAIN_nocl),file="C:\\Users\\eschumacher\\Documents\\GitHub\\ZAIN_insitu_exsitu\\Data_Files\\ZAIN_geninds\\ZAIN_garden_wild_nocl.csv")

###################################################
############ Relatedness Reduction ################
###################################################
###reduce relatedness
ZAIN_reorg_relate_df <- Demerelate(ZAIN_garden_wild_df, object = T, value = "loiselle")

##now identify how many individuals have greater than 25% relatedness = half siblings
ZAIN_halfsib_names <- names(which(unlist(ZAIN_reorg_relate_df$Empirical_Relatedness) > 0.25))

##then use this to create a document which has all the unique individual numbers for every highly related individuals
ZAIN_halfsib_names_cleanfront <- gsub("^.*\\.","", ZAIN_halfsib_names)

ZAIN_halfsib_names_cleanback <- gsub("^.*\\_","", ZAIN_halfsib_names_cleanfront)

ZAIN_relate_ind_remove <- unique(ZAIN_halfsib_names_cleanback)

##then subset genind file
butternutgen_relatedness_reduced <- butternutgen_nomd[!rownames(butternutgen_nomd@tab) %in% relate_ind_remove,]

##subset data frame 
reorg_relatedness_reduced <- reorg_relatedness[!reorg_relatedness$Ind %in% relate_ind_remove,]

###name pops
levels(butternutgen_relatedness_reduced@pop) <- butternut_24pop_names

