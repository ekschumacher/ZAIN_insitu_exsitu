##This script was created to compare nanodrop redos of ZAIN 
##There were extremely variable results in the 260/230 ratios which could indicate
##Extraction issues (high carbohyrdates remaining in the extractions)
##Or it could indicate issues with nanodropping protocols 
##This script compares 260/280 ratios and 260/230 ratios for the ZAIN individuals
##Identified to have variable ratios before and after blanks were switched out 
##To determine if the issue was extraction protocols or the blank switch 

##Libraries 

##Load in file
#set working directory
setwd("C:\\Users\\eschumacher\\Documents\\GitHub\\ZAIN_insitu_exsitu\\ZAIN_Nanodrop_Redos")
#load dataframe 
ZAIN_nano_redo <- read.csv("Data_Files\\ZAIN_nano_redo.csv")

#######################
###### Analysis #######
#######################

######Carbohydrate contamination - assessed by the 260/230 ratio
##first select data 
carb_contam <- cbind(ZAIN_nano_redo$ï..Ind, 
                     ZAIN_nano_redo %>% select(contains(c("230"))))

##then test for normality   
#plot histograms 
pdf("Figures\\histogram_260_230.pdf", width = 8, height = 6)
hist(carb_contam[,2], main = "Before Blank Switch", xlab = "260/230 Ratio")
hist(carb_contam[,3], main = "After Blank Switch", xlab = "260/230 Ratio")
dev.off()
#test normality 
shapiro.test(carb_contam[,2])
shapiro.test(carb_contam[,3])
#the "before" treatment are not normally distributed, need to use a non-parametric test
##run non-parametric paired t test 
wilcox.test(carb_contam[,2], carb_contam[,3], paired = TRUE)

##plot boxplot
pdf("Figures\\boxplot_260_230.pdf", width = 8, height = 6)
boxplot(carb_contam[,2], carb_contam[,3], ylim = c(-10, 15), 
        names = c("Before Blank Switch", "After Blank Switch"),
        xlab = "Treatment", ylab = "260/230 Ratio")
dev.off()

######RNA contamination - assessed by the 260/280 ratio 
RNA_contam <- cbind(ZAIN_nano_redo$ï..Ind, 
                    ZAIN_nano_redo %>% select(contains(c("280"))))

##then test for normality   
#plot histograms 
pdf("Figures\\histogram_260_280.pdf", width = 8, height = 6)
hist(RNA_contam[,2], main = "Before Blank Switch", xlab = "260/280 Ratio", 
     xlim = c(1,2.5))
hist(RNA_contam[,3], main = "After Blank Switch", xlab = "260/280 Ratio", 
     xlim = c(1, 2.5))
dev.off()
#test normality 
shapiro.test(RNA_contam[,2])
shapiro.test(RNA_contam[,3])
#neither data set is normally distributed (these data are left-skewed)
#need to use a non-parametric paired t-test
wilcox.test(RNA_contam[,2], RNA_contam[,3], paired = TRUE)

##Plot boxplot 
pdf("Figures\\boxplot_260_280.pdf", width = 8, height = 6)
boxplot(RNA_contam[,2], RNA_contam[,3], ylim = c(1, 2.5),
        names = c("Before Blank Switch", "After Blank Switch"),
        xlab = "Treatment", ylab = "260/280 Ratio")
dev.off()

sessionInfo()
#R version 4.0.5 (2021-03-31)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 19041)

#Matrix products: default

#locale:
#[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
#[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
#[5] LC_TIME=English_United States.1252    

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     

#loaded via a namespace (and not attached):
#[1] compiler_4.0.5 tools_4.0.5    tinytex_0.32   xfun_0.24     