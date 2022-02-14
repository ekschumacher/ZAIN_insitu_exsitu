##This is the first script for the Zamia integrifolia in situ/ex situ project. 
#We are comparing ex situ collection of wild diversity for the rare cycad Zamia 
#integrifolia. This script will focus on individual reduction based 
#on missing data and relatedness

#############################
#         Libraries         # 
#############################

library(adegenet)
library(diveRsity)
library(poppr)
library(hierfstat)
library(Demerelate)

###################################
#        Load in data files       #
###################################
#set working directory
setwd("../Data_Files")

#convert to a genepop file if necessary
#arp2gen("ZAIN_adegenet_files/ZAIN_allpop.arp")

#load in genind 
ZAIN_allpop_gen <- read.genepop("ZAIN_adegenet_files/ZAIN_allpop.gen", ncode = 3)

#ZAIN score df 
ZAIN_allpop_df <- read.csv("ZAIN_data_frames/ZAIN_allpop_df.csv")

#create a population names list 
ZAIN_allpop_names <- unique(ZAIN_allpop_df$Pop)
#name populations in the genind object  
levels(ZAIN_allpop_gen@pop) <- ZAIN_allpop_names

#name the individuals in the genind object
rownames(ZAIN_allpop_gen@tab) <- ZAIN_allpop_df$Sample.Name

########################################################
#     Removes clones and ind with missing data         #
########################################################
##Create a list of clones from the ZAIN inds 
ZAIN_clones <- as.genclone(ZAIN_allpop_gen)
##create MLG list 
ZAIN_mlg_list <- mlg.id(ZAIN_clones)
#Function to pull out individual indices where clone length greater than 1
ZAIN_clone_index <- which(sapply(ZAIN_mlg_list, function(x) length(x)>1))
##create a data frame of all of the clones 
write.csv(ZAIN_mlg_list[ZAIN_clone_index],file = "ZAIN_data_frames/ZAIN_clone_list.csv")

#This removes clones and then saves as new file for Genlex if desired
ZAIN_nocl_gen <- clonecorrect(ZAIN_allpop_gen)

##remove individuals with too much missing data 
ZAIN_allpop_clean_gen <- missingno(ZAIN_nocl_gen, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE)
##write out genepop file that's been cleaned for missing data and clones 
genind2genalex(ZAIN_allpop_clean_gen, "ZAIN_data_frames/ZAIN_allpop_clean_genalex.csv")

##write out data frame and genepop files that are cleaned for missing data and clones 
ZAIN_allpop_clean_df <- ZAIN_allpop_df[ZAIN_allpop_df$Sample.Name %in% rownames(ZAIN_allpop_clean_gen@tab),]
##write out data frame cleaned for missing data and clones 
write.csv(ZAIN_allpop_clean_df, "ZAIN_data_frames/ZAIN_allpop_clean_df.csv")

######################################
#       Relatedness Reduction        #
######################################
##write a loop to reduce data by relatedness for the garden and wild individuals separately 
ZAIN_relate_ind_remove <- list()

##list of pop type 
pop_type_list <- c("Garden", "Wild")

for(pop_type in 1:length(pop_type_list)){
  
  ##create data frames separated by type - garden or wild 
  ZAIN_df <- ZAIN_allpop_clean_df[ZAIN_allpop_clean_df$Garden_Wild == paste0(pop_type_list[[pop_type]]),]
  
  ##limit genind objects by garden of wild type 
  ZAIN_gen <- ZAIN_allpop_clean_gen[rownames(ZAIN_allpop_clean_gen@tab) %in% ZAIN_df$Sample.Name,]
  
  ##now run relatedness analysis on the data frame 
  ZAIN_relate_df <- Demerelate(ZAIN_df[,-3], object = T, value = "loiselle")
  
  ##now identify how many individuals have greater than 25% relatedness = half siblings
  ZAIN_halfsib_names <- names(which(unlist(ZAIN_relate_df$Empirical_Relatedness) > 0.25))
  
  ##create a list of related individual names 
  ZAIN_halfsib_names_cleanfront <- gsub("^.*\\.","", ZAIN_halfsib_names)
  
  ZAIN_halfsib_names_cleanback <- gsub("^.*\\_","", ZAIN_halfsib_names_cleanfront)
  
  ZAIN_relate_ind_remove[[pop_type]] <- unique(ZAIN_halfsib_names_cleanback)
  
  ##first subset the genind object with the list 
  ZAIN_rel_gen <- ZAIN_gen[!rownames(ZAIN_gen@tab) %in%  ZAIN_relate_ind_remove[[pop_type]],]
  
  #write a genalex file 
  genind2genalex(ZAIN_rel_gen, paste0("ZAIN_data_frames/ZAIN_rel_", pop_type_list[[pop_type]], "_genalex.csv"),
                 overwrite = TRUE)
  
  ##also reduce data frames by relatedness
  ZAIN_rel_df <- ZAIN_df[ZAIN_df$Sample.Name %in%  ZAIN_relate_ind_remove[[pop_type]],]
  
  ##write out data frame reduced for highly related individuals 
  write.csv(ZAIN_rel_df, paste0("ZAIN_data_frames/ZAIN_rel_", pop_type_list[[pop_type]], "_df.csv"))
  
}
