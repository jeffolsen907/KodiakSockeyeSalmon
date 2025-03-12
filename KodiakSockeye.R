#Data analysis for Kodiak Sockeye Salmon study that examines the genetic outcome of the Frazer Lake introduction
#Use base R and tidyverse functions


# Open packages--------------------------------------------------------------------------------------------------------------
library(adegenet)
library(allelematch)
library(ape)
library(cowplot)
library(genepop)
library(ggforce)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(gtools)
library(hierfstat)
library(MASS)
library(mmod)
library(pegas)
#library(PopGenKit)
library(poppr)
library(RColorBrewer)
library(reshape2)
library(spaa)
library(tidyverse)
library(viridis)
library(patchwork)
library(vegan)


#Clear workspace
rm(list = ls())

#Set print to 200 lines-----------------------------------------------------------------------------------------------------
options(tibble.print_max = 200, tibble.print_min = 200)

#Testing one, two, three


# 1=PREPARE DATA------------------------------------------------------------------------------------------------------


# 1.1-Read and arrange input file------------------------------------------------------------------------------------
#The data file "KodiakSockFrazerPaper_109.txt" was converted from genepop file 'KodiakSock42.txt' using python script.
#The new data file includes 33 of the 42 collections and recodes the collection codes using 'KodSockPopCodes.txt' (see python script Genepop2Dataframe.py).
#read input as a tibble.
input <- read_csv("KodSockFrazerPaper_109.txt") %>%
  #sort alphabaetically by Pop column.
  arrange(Pop) %>% 
  #move three mtDNA SNPs to end of tibble.
  select(-c(One_CO1, One_Cytb_17, One_Cytb_26), everything()) %>% 
  #modify mtDNA haplotype code to 3 digit from 6 digit (e.g., from 001/001 to 001).
  mutate(One_CO1 = str_sub(One_CO1, 1, 3), One_Cytb_17 = str_sub(One_Cytb_17, 1,3), One_Cytb_26 = str_sub(One_Cytb_26, 1, 3))
#modify numbers on Pop names to match numbers on map (Figure 1 of paper)
input$Pop <- plyr::mapvalues(input$Pop, unique(input$Pop), paste(str_sub(unique(input$Pop),1,3),c(1:33),sep=""))


# 1.2-Recode microsatellites-------------------------------------------------------------------------------------------
#Use "recode" in dplyr package to recode microsatellite loci from character number to allele size estimate.
#This script uses the Locus Control file 'LocusControl.txt' from ADF&G.  
LocusControl <- dget("LocusControl.txt")
#A vector of character numbers ('001', '002',....'099') that aligns with the allele size data from LocusControl.
MsatCode <- c(paste('00', c(1:9), sep = ''), paste('0', c(10:99), sep = '')) 
#the column names (Locus names) of the microsatellites in the input file.
MsatName <- names(input[96:108]) 
#Loop through each microsatellite converting each allele from character code to allele size using the mapvalue function
for (i in 1 : length(MsatName)) {
    #convert the msat list of allele sizes from list to vector.
  iMsatSize <-  unlist(LocusControl$alleles[MsatName[i]], use.names = F) 
    #add a zero to allele sizes less than 100.
  iMsatSize <- str_pad(iMsatSize, width = 3, side = 'left', pad = 0) 
    #modify the MsatCode vector so that the allele size and allele character code vectors have the same number of characters  (for mapvalues).
  iMsatCode <- MsatCode[1:length(iMsatSize)] 
    #use 'setNames' to assign MsatSize to MsatCode for recoding from allele code to allele size below
  SizeCode <- setNames(iMsatSize, iMsatCode) 
    #use 'separate' in tidyr package to create two allele columns from a single genotype column.
  input <- separate(input, !!MsatName[i], into = c('A1', 'A2'), sep = '/') 
    #use 'recode' in dplyr package to convert first allele from allele code to allele size.
  input$A1 <- recode(input$A1, !!!SizeCode) 
    #use 'recode' in dplyr package to convert second allele from allele code to allele size.
  input$A2 <- recode(input$A2, !!!SizeCode) 
    #use 'unite' in the tidyr package to merge the two recoded allele columns into a single genotype column.
  input <- unite(input, !!MsatName[i], A1, A2, sep = '/') 
}
#remove temp variables used for loop
rm("i", "iMsatCode", "iMsatSize", "LocusControl", "MsatCode", "MsatName", "SizeCode")


# 1.3-TABLE 1 : collection summary-------------------------------------------------------------------------------- 
KodiakCollections <- read_csv("KodiakCollections.csv", col_types = strrep("c", 16)) %>% 
  select(-c(Archive, USFWS_code, ADFG_code, BurgerLat, BurgerLong)) %>% 
  filter(Drainage %in% c('Karluk Lk', 'Frazer Lk', 'Ayakulik R', 'Becharof Lk')) %>% 
  filter(Code != "KBL3") %>% 
  arrange(Code) %>% 
  mutate(Code = str_sub(Code, 1,3)) %>% 
  mutate(Drainage = replace(Drainage, Drainage == "Ayakulik R", "Red")) %>% 
  mutate(Drainage = str_replace(Drainage, " Lk", "")) %>%
  mutate(Location = str_replace(Location, " Red Lk", "")) %>% 
  mutate(Habitat_Type = str_replace(Habitat_Type, "Trib_Lat", "Trib"),
         Habitat_Type = str_replace(Habitat_Type, "Trib_Term", "Trib")) %>% 
  select(-c(Frazer_Donor)) %>% 
  mutate(Map_Code = c(1:33)) %>% 
  mutate(Map_Code = replace(Map_Code, Map_Code == "13", "13d"),
         Map_Code = replace(Map_Code, Map_Code == "31", "31d"),
         Map_Code = replace(Map_Code, Map_Code == "33", "33d"))%>% 
  rename("Lake" = "Drainage", "Group" = "Code", "Timing" = "Spawn_Time", "Habitat" = "Habitat_Type", "ColNo" = "Map_Code") %>% 
  select(Lake, Location, everything())
write.table(KodiakCollections,"FrazerTable1.txt",quote=FALSE,row.names=F,col.names=T, sep = ',')


# 1.4-TABLE S2 : genotype data-----------------------------------------------------------------------------------
#Write genotype data to file using basic R.
write.table(input,"TableS2_FrazerPaper.txt",quote=FALSE,row.names=F,col.names=T,sep="\t",na="") 
#Write genotype data to file using Tidyverse readr.
write_csv(input,"TableS2_FrazerPaper.csv", col_names = T) 


# 1.5-Input subset of nuclear loci and remove samples with no score for all loci---------------------------------
#Remove mtDNA loci. Remaining 106 loci are nuclear SNPs (103) and microsatellites (13)
input_Nuc <- select(input, -c(One_CO1, One_Cytb_17, One_Cytb_26)) %>% 
  #remove rows with all '000/000' score (removes individuals with zeros at all 106 loci).
  filter_at(vars(-c(Ind,Pop)), any_vars(. != "000/000"))


# 1.6-Remove SNP loci with cumulative minor allele frequency (over all collections) less than 0.01-------------------
#define list to populate with alleles and allele frequency for each locus
AlleleList <- list()
#set minor allele frequency threshold for dropping loci
TestFreq <- 0.01
#vector of names of SNP loci 
SNPs <- names(input_Nuc[,c(3:95)])
#for loop to create summary list of tables of allele numbers and freqencies.
for (i in 1:length(SNPs)) {
  Alleles <- separate(input_Nuc, !!SNPs[i], into = c('A1', 'A2'), sep = '/')
  AlleleTable <-  table(c(Alleles$A1, Alleles$A2), exclude = '000')  
  AlleleTableProp <-  round(prop.table(AlleleTable),3)
  AlleleList[[SNPs[i]]] <-  AlleleTableProp
}
#create a subset list of SNPs from AlleleList that have minor allele frequency less than TestFreq (e.g., 0.01) 
DropSNPs <- subset(AlleleList, lapply(AlleleList, function(x) min(x) < TestFreq) == TRUE) 
#remove SNPs with cummulation MAF less than 0.01 
input_Nuc <- select(input_Nuc, -c(names(DropSNPs)))
#remove temp variables used for loop
rm("i", "AlleleList", "Alleles", "AlleleTable", "AlleleTableProp", "TestFreq", "SNPs")


# 1.7-Create input file for mtDNA SNPs, remove samples with no score and recode to single haplotype------------------------
  #use 'select' in dplyr package to create tibble with only mtDNA loci.
input_mtDNA <- select(input, c(Ind, Pop, One_CO1, One_Cytb_17, One_Cytb_26)) %>% 
  #remove rows with all '000' score (removes individuals with zeros at all 3 loci).
  filter_at(vars(One_CO1, One_Cytb_17, One_Cytb_26), all_vars(. != "000")) %>% 
  #add column of single concatenated haplotype
  mutate(HaploCode = paste(str_sub(.$One_CO1,3), str_sub(.$One_Cytb_17,3), str_sub(.$One_Cytb_26,3), sep = "")) %>% 
  #remove individual haplotypes and retain combined haplotype 
  select(c(Ind, Pop, HaploCode))


# 1.8-Prepare nuclear and mtDNA data for use in adegenet, hierfsat and poppr packages--------------------------------------------

#Nuclear loci  
#Vector of populations from nuclear data (for adegenet)
PopList <- input_Nuc$Pop 
#Vector of individuals from nuclear data (for adegenet)
IndList <- input_Nuc$Ind 
#Vector of loci from nuclear data
LociList <- names(select(input_Nuc, -c(Pop, Ind)))
#Use 'select' in dplyr package to remove columns Ind, Pop.
input4AD <- select(input_Nuc, -c(Ind, Pop)) 
#Convert to adegent input for nuclear loci
AD <- df2genind(input4AD,sep = '/', ncode=6, ind.names=IndList, pop=PopList, NA.char="000/000", type="codom") 
#Convert to hierfsat input for nuclear loci
HF <- genind2hierfstat(AD) 

#mtDNA
#Vector of populations for Adegent
PopList_mt <- input_mtDNA$Pop 
#Vector of individuals from Adegent
IndList_mt <- input_mtDNA$Ind 
#Vector of mtDNA haplotypes
LociList_mt <- names(select(input_mtDNA, -c(Pop, Ind)))
#Use 'select' in dplyr package to remove columns Ind, Pop.
input4AD_mtDNA <- select(input_mtDNA, -c(Ind, Pop))
#Convert to adegent input for mtDNA haplotypes
AD_mtDNA = df2genind(input4AD_mtDNA,sep = NULL, ncode=3, ind.names = IndList_mt, pop = PopList_mt,NA.char="000",type="codom")
#Convert to hierfsat input for mtDNA haplotypes
HF_mtDNA = genind2hierfstat(AD_mtDNA) 


#1.9-Prepare nuclear data for use in genepop package--------------------------------------------------------
#remove "One_' from prefix of each locus name (so can distinguish loci in output files from GenePop analysis)
names(input4AD) <- gsub("One_","",names(input4AD),fixed=T) 
#Creates a dataframe with Individual and Population in first two columns followed by columns for each locus in three digit format ('010/030').  
inputGP <-  cbind(IndList, PopList, input4AD) 
#Creates blank dataframe with same number of columns as inputGP and no. columns - 1 rows.
NewGP <-  data.frame(matrix(nrow = length(inputGP) - 1, ncol = length(inputGP))) 
#Assign the first column of NewGP the column names of inputGP column 2 through end.  
NewGP[,1] <-  names(inputGP[,2:length(inputGP)]) 
#Assign the column names of inputGP to NewGP
names(NewGP) <-  names(inputGP) 
#Assign the first row, column 1 of NewGP the string "Genepop Input'.
NewGP[1,1] <- 'Genepop Input' 
#Create LociCount vector with one element that is the number of columns of NewGP minus 2.
LociCount <-  length(NewGP)-2 
#for loop replaces/removes all slashes separating alleles at each locus.
for(i in 1:LociCount) {
  inputGP[,i+2] = gsub('/', '', inputGP[,i+2])
}
#Add a comma to each character in the column IndList in the dataframe inputGP.  The comma in front of the individual ID is needed for genepop. 
inputGP$IndList <-  paste(inputGP$IndList, ',', sep = '') 
#Create a one row empty dataframe 
Pop <-  as.data.frame(t(c(rep(NA, length(inputGP))))) 
#Assign the colum names of inputGP to Pop
names(Pop) <-  names(inputGP) 
#Assign the first row, column 1 of Pop the string 'Pop'
Pop[1] <-  'Pop' 
#Create a PopN vector with the names of each population from the the column inputGP$PopList in inputGP.
PopN <-  unique(inputGP$PopList) 
#for loop appends the datafram NewGP with the Pop dataframe and each population, in turn, from the dataframe inputGP.
for(i in 1:length(PopN)) {
  tempGP <-  subset(inputGP, inputGP$PopList == PopN[i])
  NewGP <-  rbind(NewGP, Pop, tempGP)  
}
#Remove the column titled PopList from NewGP
NewGP <-  subset(NewGP, select = -c(PopList)) 
#Write NewGP to file "KodiakGPFrazer.txt.
write.table(NewGP,"KodiakGPFrazer.txt",quote=FALSE,row.names=F,col.names=F,sep="\t",na="") 
#remove temp variables used for loop
rm("i", "inputGP", "LociCount", "NewGP", "Pop", "PopN", "tempGP")



# 2=FIRST LEVEL ANALYSIS----------------------------------------------------------------------------------------------


# 2.1-hierFstat : Ho, Hs, Ar, Fis, Fst----------------------------------------------------------------------------------
#For all nuclear loci and populations
B <-  basic.stats(HF)
#compute Hs, Fst for all mtDNA and populations
B_mtDNA = basic.stats(HF_mtDNA, diploid = F) #use hierFstat to compute Hs for all mtDNA and populations
#Compute Ar for all nuclear loci and populations
AR <- allelic.richness(HF) 
#Use "select" in the dplyr package to create SNP subset from hierFstat input file for nuclear loci 
HF_SNPs <- select(HF, pop, One_ACBP_79:One_zP3b) 
#Compute Ho, Hs, Fis, Fst for all SNP loci and populations
B_SNPs <- basic.stats(HF_SNPs) 
#use "select" in the dplyr package to create mSat subset from hierFstat input file for nuclear loci 
HF_mSats <- select(HF, pop, One_Omy77v4:One_uSat60v4) 
#use hierFstat to compute Ho, Hs, Fis, Fst for all mSat loci and populations
B_mSats <- basic.stats(HF_mSats) 
#Compute Weir and Cockerham estimates of Fst for all SNP loci and populations
WC_SNPs <- wc(HF_SNPs) 
#Compute Weir and Cockerham estimates of Fst for all mSat loci and populations
WC_mSats <- wc(HF_mSats) 
#Compute Ar for all SNP loci and populations
AR_SNPs <- allelic.richness(HF_SNPs) 
#Compute Ar for all mSat loci and populations
AR_mSats <- allelic.richness(HF_mSats) 
#compute mean Ar for each SNP locus across populations
AR_SNPs$mean <- round(rowMeans(AR_SNPs$Ar),2) 
#compute mean Ar for each mSat locus across populations
AR_mSats$mean <- round(rowMeans(AR_mSats$Ar),2) 


# 2.2-Poppr : Shannon-Weiner index---------------------------------------------------------------------
# for each nuclear locus - used to select most informative locus of locus pair if linked
SW_Nuc <- locus_table(AD, index = 'shannon')  
#Coerce SW_Nuc to tibble
SW_Nuc <- as_tibble(rownames_to_column(as.data.frame(SW_Nuc), var = "Locus")) 
#Use "filter" in dplyr package to remove last row (says "mean" in Locus column) showing mean SW estimates for all loci
SW_Nuc <- filter(SW_Nuc, Locus != "mean") 


# 2.3-TABLE S1 : descriptive statistics for all loci--------------------------------------------
#Create dataframe of descriptive stats for all loci from hierFstat and poppr output above
ByLoc <-  rbind(B$perloc, B_mtDNA$perloc[1,], B_SNPs$overall, B_mSats$overall) %>% 
  rownames_to_column(var = "Locus") %>% 
  #coerce ByLoc to tibble
  as.tibble() %>% 
  #create columns of statistics for each locus
  mutate(Class = c(rep("SNP",90), rep("mSat",13), "SNPmt", "allSNP", "allmSat"), 
         Ref = c(rep("a", 106)),
         FstWC = c(WC_SNPs$per.loc$FST, WC_mSats$per.loc$FST, B_mtDNA$perloc[1,c("Fst")], WC_SNPs$FST, WC_mSats$FST),
         FisWC = c(WC_SNPs$per.loc$FIS, WC_mSats$per.loc$FIS, B_mtDNA$perloc[1,c("Fis")], WC_SNPs$FIS, WC_mSats$FIS),
         AR = c(AR_SNPs$mean, AR_mSats$mean, NA , mean(AR_SNPs$mean), mean(AR_mSats$mean)),
         SW = c(SW_Nuc$H, NA, NA, NA)) 
#use "select" in dplyr package to create subset table of summary statistics for each locus
ByLocSub <- select(ByLoc, Locus, Class, Ref, Ho, Hs, FisWC, FstWC, AR, SW) 
#use "mutate_if" in dplyr package to round number to 3 decimals
ByLocSub <- mutate_if(ByLocSub, is.double, round, 3) 
#write "ByLocSub" to file using Tidyverse readr
write_csv(ByLocSub, "TableS1_FrazerPaper.csv") 


# 2.4-TABLE S3 : descriptive statistics for all collections---------------------------------------
#create dataframe of descriptive stats from hierFstat and poppr output above
#create vector of mean Hs for mSats for each population
mSat_Hs <- round(colMeans(B_mSats$Hs),3) 
#create vector of mean Hs for SNPs for each population
SNP_Hs <- round(colMeans(B_SNPs$Hs),3) 
#create vector of mean Hs for all nuclear loci for each population
Allnuc_Hs <- round(colMeans(B$Hs),3) 
#create vector of mean Ar for mSats for each population
mSat_Ar <- round(colMeans(AR_mSats$Ar),1) 
#create vector of mean Ar for SNPs for each population
SNP_Ar <- round(colMeans(AR_SNPs$Ar),1) 
#create vector of mean Ar for all nuclear loci for each population
Allnuc_Ar <- round(colMeans(AR$Ar),1) 
#create vector of Hs for mtDNA for each population
mtDNA_Hs <- round(B_mtDNA$Hs[1,],3) 
#create dataframe from vectors above of diversity statistics for each population
PopsDiversity <- cbind(mSat_Hs, SNP_Hs, Allnuc_Hs, mtDNA_Hs, mSat_Ar, SNP_Ar, Allnuc_Ar) 
#coerce "PopsDiversity" to tibble
PopsDiversityTib <- as.tibble(rownames_to_column(as.data.frame(PopsDiversity), var = "Pop")) 
#use 'mutate' to add column with Lake name for each collection
PopsDiversityTib <- mutate(PopsDiversityTib, Lake = c(rep("Frazer", 10), rep("Karluk", 16), rep("Red", 6), "Ruth"))
#use 'select' to move column 'Lake' to beginning
PopsDiversityTib <- select(PopsDiversityTib, Lake, everything())
#use aggregate to compute mean values for Hs, Ar for each Lake
aggregate(select(PopsDiversityTib, mSat_Hs:Allnuc_Ar), by = list(PopsDiversityTib$Lake), FUN = mean)
#Write "PopsDiversityTib" to file using Tidyverse readr
write_csv(PopsDiversityTib,"TableS3_FrazerPaper.csv", col_names = T)

#Compute max and min Hs for SNPs and mSats
max(B_SNPs$Hs) #max Ho for SNPs
min(B_SNPs$Hs) #min Ho for SNPs
max(B_mSats$Hs)  #max Ho for mSats
min(B_mSats$Hs)  #min Ho for mSats


# 2.5-Genepop : HWP and GD----------------------------------------
#Define locinfile from geneopop input file for nuclear loci created above
locinfile <-  "KodiakGPFrazer.txt"
#Compute allele and genotype frequencies for nuclear loci
basic_info(locinfile, outputFile = "KodiakGPFrazerOut.txt")
#Test HWP (Hardy_Weinberg proportions) for all nuclear loci and collections
test_HW(locinfile, outputFile = 'KodiakGPFrazerHW.txt', enumeration = TRUE, dememorization = 10000, batches = 100, iterations = 5000)
#Test GD (gametic disequilibrium) for all pairs of nuclear loci across collections
test_LD(locinfile, outputFile = 'KodiakGPFrazerGD.txt', dememorization = 10000, batches = 100, iterations = 5000)

# 2.6-Evaluate results of HWP tests from genepop------------------------------
#First, test assumption of global HWP across all loci and collections
#use cumulative binomial probability distribution to derive 95% CI to the number of table-wide significant tests (p < 0.05). See Waples 2015. 
#Modify the test_HW output using python script (outside of R) to create a locus x collection matrix of p-values. 
#Read data file matrix (created using python script) of HWP p-values 33 pops. The col_types argument explicitly specifies the column data type (c=character, d=double).
inputHWP <- read_csv('FrazerGenepopHWP_Pvalues.txt', col_types = paste('c',strrep('d',33), sep=''), na = '-') 
#Create tibble from inputHWP that sums the number p-values < 0.05 for each locus.
LocLTalfa <- tibble('Loc' = inputHWP$Loc, 'NLocLTalfa' = rowSums(inputHWP[-1] < 0.05, na.rm = TRUE)) 
#Create tibble from inputHWP that sums the number p-values < 0.05 for each population.
PopLTalfa <- tibble('Pop' = names(inputHWP[-1]), 'NPopLTalfa' = colSums(inputHWP[-1] < 0.05, na.rm=TRUE)) 
#table-wide number of HWP tests (number of loci x number of populations)
NoTestsHWP <- length(LocLTalfa$Loc)*length(PopLTalfa$Pop) 
NoTestsSigHWP <- sum(LocLTalfa$NLocLTalfa)
binomtestHWP <- dbinom(1:200, size = NoTestsHWP, prob = 0.05)
binomtestHWP <- round(binomtestHWP,4)
  #cumulative binomial probability distribution to determine if the number of significant tests N (e.g., N = sum(LocLTalfa$NLocLTalfa)) is within (NS) the 0.025 < N < 0.975 interval.
binomtestcumHWP <- cumsum(binomtestHWP) 

#Second, test assumption of HWP for each collection across loci and each locus across collections

#Population-level test of HWP across loci
PopSigTest <-  read_csv("FrazerGenepopHWP_PopSigTest.txt", col_types = paste('c',strrep('d',6), sep=''), na = '-')
PopBinomTest <- dbinom(0:length(LocLTalfa$Loc), size = length(LocLTalfa$Loc), prob = 0.05)
PopBinomTestTable <- tibble('NumSigLoci' = 0:length(LocLTalfa$Loc), 'PopBinomExp' = round(PopBinomTest * length(PopLTalfa$Pop), 2), 'PopBinomObs' = as.vector(table(factor(PopSigTest$ST, levels = 0:length(LocLTalfa$Loc)))))
PopBinomTestTableTidy <- gather(PopBinomTestTable, 'PopBinomExp', 'PopBinomObs', key = 'ObsOrExp', value = 'NumPop')
#Plot results for collections across loci
PopBinomTestPlot <- ggplot(PopBinomTestTableTidy) + 
  geom_bar(aes(NumSigLoci, NumPop, fill = ObsOrExp), position = 'dodge',  stat = 'identity') + 
  scale_fill_manual(values = c("black", "gray"), labels = c(" Expected", " Observed")) + 
  theme_bw() + 
  ggtitle("Collections") +
  scale_y_continuous(expand = c(0,0), limits = c(0,8)) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0,103,10)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.2, 0.8), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.margin=unit(c(0.5,10,0.5,1.2),"cm"),
        legend.title = element_blank(), legend.key.size = unit(0.3,"cm"),
        plot.title = element_text(vjust = -10, hjust = 0.5)) +
  annotate("text", x = 21, y = 1.5, label = "KTE18", size = 3)

#Locus-level test of HWP across populations
LocSigTest <-  read_csv("FrazerGenepopHWP_LocSigTest.txt", col_types = paste('c',strrep('d',6), sep=''), na = '-')
LocBinomTest <- dbinom(0:length(PopLTalfa$Pop), size = length(PopLTalfa$Pop), prob = 0.05)
LocBinomTestTable <- tibble('NumSigCollections' = 0:length(PopLTalfa$Pop), 'LocBinomExp' = round(LocBinomTest * length(LocLTalfa$Loc), 2), 'LocBinomObs' = as.vector(table(factor(LocSigTest$ST, levels = 0:length(PopLTalfa$Pop)))))
LocBinomTestTableTidy <- gather(LocBinomTestTable, 'LocBinomExp', 'LocBinomObs', key = 'ObsOrExp', value = 'NumLoci')
#Plot results for loci across collections
LocBinomTestPlot <- ggplot(LocBinomTestTableTidy) + 
  geom_bar(aes(NumSigCollections, NumLoci, fill = ObsOrExp), position = 'dodge',  stat = 'identity') + 
  scale_fill_manual(values = c("black", "gray"), labels = c(" Expected", " Observed")) + 
  theme_bw() + 
  ggtitle("Loci") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,40,5), limits = c(0,40)) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0,33,5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.2,0.8), 
        axis.title = element_text(size=12), axis.text = element_text(size=10), plot.margin=unit(c(0.5,10,0.5,1.2),"cm"),
        legend.title = element_blank(), legend.key.size = unit(0.3,"cm"),
        plot.title = element_text(vjust = -10, hjust = 0.5)) +
  annotate("text", x = 11, y = 3, label = "One_ACBP_79", size = 3)


# 2.7-FIGURE S1 : Plot binomial distribution of HWP results for collections and loci---------------------------------------------- 
tiff('FigureS1_FrazerPaper.tiff', height = 6.5, width = 8.5, units = 'in', compression = 'none', res = 150) 
ggarrange(PopBinomTestPlot, LocBinomTestPlot, ncol=1, nrow=2, align = 'v')
dev.off()


# 2.8-Evaluate results of GD tests from genepop-----------------------------------------------------------------------------------
#read data file matrix of HWP p-values 33 pops
inputGD <- read_csv("FrazerGenepopGD_Pvalues.txt", col_types = paste('c', strrep('d',33), strrep('i',2), sep=''), na='NA') 
NoTestsGD <- length(inputGD$LocXLoc)*length(PopLTalfa$Pop)
NoTestsSigGD <- sum(inputGD$ST)
binomtestGD <- dbinom(1:10000, size = NoTestsGD, prob = 0.05)
binomtestGD <- round(binomtestGD,4)
#cumulative binomial probability distribution to determine if the number of significant tests N (e.g., N = sum(LocLTalfa$NLocLTalfa)) is within (NS) the 0.025 < N < 0.975 interval.
binomtestcumGD <- cumsum(binomtestGD) 
arrange(select(inputGD, LocXLoc, ST), desc(ST))

GDBinomTest <- dbinom(0:length(PopLTalfa$Pop), size = length(PopLTalfa$Pop), prob = 0.05)
GDBinomTestTable <- tibble('NumSigCollections' = 0:length(PopLTalfa$Pop), " Expected" = round(GDBinomTest * length(inputGD$LocXLoc), 2), " Observed" = as.vector(table(factor(inputGD$ST, levels = 0:length(PopLTalfa$Pop)))))
GDBinomTestTableTidy <- gather(GDBinomTestTable, " Expected", " Observed", key = 'ObsOrExp', value = 'NumLocPairs')

#plot of count of locus pairs versus number of significant collections
GDBinomTestPlot <- ggplot(GDBinomTestTableTidy) + 
  geom_bar(aes(NumSigCollections, NumLocPairs, fill = ObsOrExp), position = 'dodge',  stat = 'identity') + 
  scale_fill_manual(values = c("black", "gray")) +
  theme_bw() + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0,33,5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        axis.title = element_blank(), axis.text = element_text(size=8), plot.margin=unit(c(0.5,10,0.5,1.2),"cm"),
        legend.title = element_blank(), plot.background = element_rect(fill = "transparent",colour = NA))
#Zoom in GD binomial plot to show locus pairs with GD in many collections 
GDBinomTestPlotZoom <- GDBinomTestPlot + coord_cartesian(ylim = c(0,40)) + 
  ggtitle("Locus Pairs") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.3,0.8),
        axis.title = element_text(size=14), axis.text = element_text(size=10), plot.margin=unit(c(0.5,10,0.5,1.2),"cm"),
        legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=10), 
        plot.title = element_text(size = 14, vjust = -10, hjust = 0.5)) +
  annotate("text", x = c(24,30), y = 3, label = c("One_MHC2_190 x\nOne_MHC2_251", "One_Tf_ex10_750 x\nOne_Tf_ex3_182"), size = 4)
#Inset overall plot within zoom plot
GDPlotInsetZoom <- ggdraw() +
  draw_plot(GDBinomTestPlotZoom) +
  draw_plot(GDBinomTestPlot, x = 0.45, y = 0.45, width = 0.5, height = 0.4)


# 2.9-FIGURE S2 : Plot of GD test results - number of locus pairs versus number of significant collections-------------------------
tiff('FigureS2_FrazerPaper.tiff', height = 8.5, width = 14, units = 'in', compression = 'none', res = 100) 
print(GDPlotInsetZoom)
dev.off()


# 3=UPDATE DATA-------------------------------------------------------------------------------------------------------

# 3.1-Remove collection KTE3 and locus One_ACBP_79 based on results of HWP tests and locus One_Tf_ex10_750 based on results of GD tests---------------
# 3.1.1-Nuclear loci--------------------------------------
input_NucGo <- select(input_Nuc, -c(One_ACBP_79, One_Tf_ex10_750)) %>% 
  filter(Pop != "KTE18")
# 3.1.2-mtDNA------------------------------------------------
input_mtDNA_Go <- filter(input_mtDNA, Pop != "KTE18")
input_mtDNA_Go4Fstat <- mutate(input_mtDNA_Go, HaploCode2 = paste(input_mtDNA_Go$HaploCode,input_mtDNA_Go$HaploCode, sep = "")) %>% 
  select(-HaploCode)
input_mtDNA_GoHapSep <- mutate(input_mtDNA_Go, HaploCode2 = paste(input_mtDNA_Go$HaploCode,input_mtDNA_Go$HaploCode, sep = "/")) %>% 
  select(-HaploCode)
# 3.1.3-nuclear loci and mtDNA-----------------------------------------
#Use "left_join" to join the mtDNA input to the nuclear loci input and retain all rows from nuclear data set.
input_NucMtDNA_GoL <- left_join(input_NucGo,input_mtDNA_GoHapSep) %>% 
  replace_na(list(HaploCode2 = "000/000"))
#Save updated input file to project folder
saveRDS(input_NucMtDNA_GoL, "input_NucMtDNA_GoL.rds")


input_NucMtDNA_GoL <- readRDS("input_NucMtDNA_GoL.rds")
input_NucGo <- select(input_NucMtDNA_GoL, -c(HaploCode2))


# 3.2-Prepare nuclear and mtDNA data for use in adegenet, hierfsat packages---------------------------------------------------------

#Nuclear loci  
#Vector of populations from nuclear data (for adegenet)
PopListNucGo <- input_NucGo$Pop #List of populations for Adegent
#Vector of individuals from nuclear data (for adegenet)
IndListNucGo <- input_NucGo$Ind #List of individuals from Adegent
#Vector of nuclear loci
LociListNucGo <- names(select(input_NucGo, -c(Pop, Ind)))
#Use 'select' in dplyr package to remove columns Ind, Pop.
inputNucGo4AD <- select(input_NucGo, -c(Ind, Pop)) 
#Convert to adegenet input file
ADnucGo <- df2genind(inputNucGo4AD,sep = '/', ncode=6, ind.names=IndListNucGo, pop=PopListNucGo, NA.char="000/000", type="codom") 
#Convert to hierfstat input file
HFnucGo <- genind2hierfstat(ADnucGo)

#mtDNA
#Vector of populations from mtDNA data (for adegenet)
PopList_mtGo <- input_mtDNA_Go$Pop 
#Vector of individuals from mtDNA data (for adegenet)
IndList_mtGo <- input_mtDNA_Go$Ind 
#Vector of mtDNA haplotypes
LociList_mtGo <- names(select(input_mtDNA_Go, -c(Pop, Ind)))
#Use 'select' in dplyr package to remove columns Ind, Pop.
input4AD_mtDNAGo <- select(input_mtDNA_Go, -c(Ind, Pop))
#Convert to adegenet input file
AD_mtDNAGo = df2genind(input4AD_mtDNAGo,sep = NULL, ncode=3, ploidy = 1, ind.names = IndList_mtGo, pop = PopList_mtGo,NA.char="000",type="codom")
#Convert to hierfstat input file
HF_mtDNAGo = genind2hierfstat(AD_mtDNAGo)


# 3.3-Separate SNP and microsatellite input files for Red, Karluk and Frazer Lakes.  Needed for outlier loci analysis--------------------------------
RedKarluk <- unique(PopListNucGo)[c(11:31)]
Frazer <- unique(PopListNucGo)[c(1:10)]
Red <- unique(PopListNucGo)[c(26:31)]
Karluk <- unique(PopListNucGo)[c(11:25)]
input_NucGoSNP <- select(input_NucGo, Ind, Pop, One_agt_132:One_zP3b)
input_NucGoSNP_RK <- filter(input_NucGoSNP, Pop %in% RedKarluk)
input_NucGoSNP_F <- filter(input_NucGoSNP, Pop %in% Frazer)
input_NucGoSNP_R <- filter(input_NucGoSNP, Pop %in% Red)
input_NucGoSNP_K <- filter(input_NucGoSNP, Pop %in% Karluk)
input_NucGoMsat <- select(input_NucGo, Ind, Pop, One_Omy77v4:One_uSat60v4)
input_NucGoMsat_RK <- filter(input_NucGoMsat, Pop %in% RedKarluk)
input_NucGoMsat_F <- filter(input_NucGoMsat, Pop %in% Frazer)
input_NucGoMsat_R <- filter(input_NucGoMsat, Pop %in% Red)
input_NucGoMsat_K <- filter(input_NucGoMsat, Pop %in% Karluk)

#SNPs
PopListNucGoSNP <- input_NucGoSNP$Pop 
IndListNucGoSNP <- input_NucGoSNP$Ind 
LociListNucGoSNP <- names(select(input_NucGoSNP, -c(Pop, Ind)))
inputNucGoSNP_4AD <- select(input_NucGoSNP, -c(Ind, Pop)) 
ADnucGoSNP <- df2genind(inputNucGoSNP_4AD,sep = '/', ncode=6, ind.names=IndListNucGoSNP, pop=PopListNucGoSNP, NA.char="000/000", type="codom") 
HFnucGoSNP <- genind2hierfstat(ADnucGoSNP) 

#Microsatellites
PopListNucGoMsat <- input_NucGoMsat$Pop 
IndListNucGoMsat <- input_NucGoMsat$Ind 
LociListNucGoMsat <- names(select(input_NucGoMsat, -c(Pop, Ind)))
inputNucGoMsat_4AD <- select(input_NucGoMsat, -c(Ind, Pop)) 
ADnucGoMsat <- df2genind(inputNucGoMsat_4AD,sep = '/', ncode=6, ind.names=IndListNucGoMsat, pop=PopListNucGoMsat, NA.char="000/000", type="codom")
HFnucGoMsat <- genind2hierfstat(ADnucGoMsat) 

#SNPs Red and Karluk
PopListNucGoSNP_RK <- input_NucGoSNP_RK$Pop 
IndListNucGoSNP_RK <- input_NucGoSNP_RK$Ind 
LociListNucGoSNP_RK <- names(select(input_NucGoSNP_RK, -c(Pop, Ind)))
inputNucGoSNP_RK4AD <- select(input_NucGoSNP_RK, -c(Ind, Pop)) 
ADnucGoSNP_RK <- df2genind(inputNucGoSNP_RK4AD,sep = '/', ncode=6, ind.names=IndListNucGoSNP_RK, pop=PopListNucGoSNP_RK, NA.char="000/000", type="codom") 
HFnucGoSNP_RK <- genind2hierfstat(ADnucGoSNP_RK) 

#SNPs Frazer
PopListNucGoSNP_F <- input_NucGoSNP_F$Pop 
IndListNucGoSNP_F <- input_NucGoSNP_F$Ind 
LociListNucGoSNP_F <- names(select(input_NucGoSNP_F, -c(Pop, Ind)))
inputNucGoSNP_F4AD <- select(input_NucGoSNP_F, -c(Ind, Pop)) 
ADnucGoSNP_F <- df2genind(inputNucGoSNP_F4AD,sep = '/', ncode=6, ind.names=IndListNucGoSNP_F, pop=PopListNucGoSNP_F, NA.char="000/000", type="codom") 
HFnucGoSNP_F <- genind2hierfstat(ADnucGoSNP_F) 
saveRDS(LociListNucGoSNP_F, "LociListNucGoSNP_F.rds")

#SNPs Karluk
PopListNucGoSNP_K <- input_NucGoSNP_K$Pop 
IndListNucGoSNP_K <- input_NucGoSNP_K$Ind 
LociListNucGoSNP_K <- names(select(input_NucGoSNP_K, -c(Pop, Ind)))
inputNucGoSNP_K4AD <- select(input_NucGoSNP_K, -c(Ind, Pop)) 
ADnucGoSNP_K <- df2genind(inputNucGoSNP_K4AD,sep = '/', ncode=6, ind.names=IndListNucGoSNP_K, pop=PopListNucGoSNP_K, NA.char="000/000", type="codom")
HFnucGoSNP_K <- genind2hierfstat(ADnucGoSNP_K) 
saveRDS(LociListNucGoSNP_K, "LociListNucGoSNP_K.rds")

#SNPs Red
PopListNucGoSNP_R <- input_NucGoSNP_R$Pop 
IndListNucGoSNP_R <- input_NucGoSNP_R$Ind 
LociListNucGoSNP_R <- names(select(input_NucGoSNP_R, -c(Pop, Ind)))
inputNucGoSNP_R4AD <- select(input_NucGoSNP_R, -c(Ind, Pop)) 
ADnucGoSNP_R <- df2genind(inputNucGoSNP_R4AD,sep = '/', ncode=6, ind.names=IndListNucGoSNP_R, pop=PopListNucGoSNP_R, NA.char="000/000", type="codom")
HFnucGoSNP_R <- genind2hierfstat(ADnucGoSNP_R) 
saveRDS(LociListNucGoSNP_R, "LociListNucGoSNP_R.rds")

#Microsatellites Frazer
PopListNucGoMsat_F <- input_NucGoMsat_F$Pop
IndListNucGoMsat_F <- input_NucGoMsat_F$Ind 
LociListNucGoMsat_F <- names(select(input_NucGoMsat_F, -c(Pop, Ind)))
inputNucGoMsat_F4AD <- select(input_NucGoMsat_F, -c(Ind, Pop)) 
ADnucGoMsat_F <- df2genind(inputNucGoMsat_F4AD,sep = '/', ncode=6, ind.names=IndListNucGoMsat_F, pop=PopListNucGoMsat_F, NA.char="000/000", type="codom") 
HFnucGoMsat_F <- genind2hierfstat(ADnucGoMsat_F) 


# 3.4-Use Arlequin outside R to test for Outlier loci-----------------------------------
#The six SNP and microsatellite data sets for Frazer, Karluk and Red lakes (above) were converted to genepop files using code above.
#These six genepop files were converted to arlequin file using arlequin v3.5.
#Each data set was tested for outlier loci using arlequin v3.5
#The arlequin results for SNPs from the three lakes were read below and plotted.
#The microsatellite results were not plotted because no outlier loci were detected.
#Read the SNP locus list for each lake (saved as rds above)
LociListNucGoSNP_K <- readRDS("LociListNucGoSNP_K.rds")
LociListNucGoSNP_R <- readRDS("LociListNucGoSNP_R.rds")
LociListNucGoSNP_F <- readRDS("LociListNucGoSNP_F.rds")

# 3.4.1-SNP outliers for Karluk Lake------------------------------------------------
SNP_K_basic <- basic.stats(HFnucGoSNP_K)
#check max and min Hs across all loci and collections (for use in Arlequin outlier loci simulation)
max(SNP_K_basic$Hs)
min(SNP_K_basic$Hs)
#use output from Arlequin
OLinputK <- read_tsv("fdist2_ObsOut_K.txt") %>% 
  select(c(1:5)) %>% 
  rename_at(vars(colnames(.)), ~ c("LocusNo", "ObsH_BP", "ObsFst", "FstPval", "One_FstQuant")) %>% 
  bind_cols(., as_tibble(gsub("One_","",LociListNucGoSNP_K,fixed=T))) %>% 
  rename(LocusName = value) %>% 
  mutate(Lake = c("Karluk")) %>% 
  select(Lake, LocusNo, LocusName, everything())
OLSimInputK <- read_table("FstvHoNull_SNP_K.txt") %>% 
  mutate(Lake = c("Karluk")) %>% 
  select(Lake, everything())

# 3.4.2-SNP outliers for Red Lake--------------------------------------------------------
SNP_R_basic <- basic.stats(HFnucGoSNP_R)
#check max and min Hs across all loci and collections (for use in Arlequin outlier loci simulation)
max(SNP_R_basic$Hs)
min(SNP_R_basic$Hs)
#use output from Arlequin
OLinputR <- read_tsv("fdist2_ObsOut_R.txt") %>% 
  select(c(1:5)) %>% 
  rename_at(vars(colnames(.)), ~ c("LocusNo", "ObsH_BP", "ObsFst", "FstPval", "One_FstQuant")) %>% 
  bind_cols(., as_tibble(gsub("One_","",LociListNucGoSNP_R,fixed=T))) %>% 
  rename(LocusName = value) %>% 
  mutate(Lake = c("Red")) %>% 
  select(Lake, LocusNo, LocusName, everything())  
OLSimInputR <- read_table("FstvHoNull_SNP_R.txt") %>% 
  mutate(Lake = c("Red")) %>% 
  select(Lake, everything())  

# 3.4.3-SNP outliers for Frazer Lake---------------------------------------------------------
SNP_F_basic <- basic.stats(HFnucGoSNP_F)
#check max and min Hs across all loci and collections (for use in Arlequin outlier loci simulation)
max(SNP_F_basic$Hs)
min(SNP_F_basic$Hs)
#use output from Arlequin
OLinputF <- read_tsv("fdist2_ObsOut_F.txt") %>% 
  select(c(1:5)) %>% 
  rename_at(vars(colnames(.)), ~ c("LocusNo", "ObsH_BP", "ObsFst", "FstPval", "One_FstQuant")) %>% 
  bind_cols(., as_tibble(gsub("One_","",LociListNucGoSNP_F,fixed=T))) %>% 
  rename(LocusName = value) %>% 
  mutate(Lake = c("Frazer")) %>% 
  select(Lake, LocusNo, LocusName, everything())  
OLSimInputF <- read_table("FstvHoNull_SNP_F.txt") %>% 
  mutate(Lake = c("Frazer")) %>% 
  select(Lake, everything())  

# 3.4.4-single tibble combining outlier loci results from each lake------------------------------------
OLinputKRF <- bind_rows(OLinputK, OLinputR, OLinputF) %>% 
  mutate(Lake = paste(Lake, "Lake"))
OLinputKRF$Lake <- factor(OLinputKRF$Lake, levels = c("Frazer Lake", "Karluk Lake", "Red Lake"))
OLSimInputKRF <- bind_rows(OLSimInputK, OLSimInputR, OLSimInputF) %>% 
  mutate(Lake = paste(Lake, "Lake"))
OLSimInputKRF$Lake <- factor(OLSimInputKRF$Lake, levels = c("Frazer Lake", "Karluk Lake", "Red Lake"))


# 3.4.5-FIGURE 2 : plot outlier loci results for SNPs for each lake------------------------------------
#creates a high resolution TIFF file for publication
Figure2 <- ggplot(OLinputKRF) + 
  geom_point(data = filter(OLinputKRF, FstPval > 0.01), aes(ObsH_BP, ObsFst), shape = 1) + 
  geom_point(data = filter(OLinputKRF, FstPval < 0.01 & ObsH_BP > 0), aes(ObsH_BP, ObsFst), shape = 2, size = 2) +
  theme_bw() +
  scale_y_continuous(name = "Fst", breaks = seq(0,0.6,0.1)) + 
  scale_x_continuous(name = "Ho/(1-Fst)", breaks = seq(0,0.7,0.1)) +
  geom_line(aes(Het_BP, `0.01`), linetype = "dashed", data = OLSimInputKRF) + 
  geom_line(aes(Het_BP, `0.99`), linetype = "dashed", data = OLSimInputKRF) + 
  coord_cartesian(xlim = c(0, 0.7), ylim = c(-0.05,0.60)) +
  geom_text(data = filter(OLinputKRF, FstPval < 0.01 & ObsH_BP > 0), 
            aes(ObsH_BP, ObsFst, label = "One_U1004_183"), nudge_x = -0.075, nudge_y = 0.02, size = 3) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin=unit(c(0,0,0,0.5),"cm"), 
        axis.title = element_text(size=10), axis.text = element_text(size=10), 
        plot.title = element_text(size=16, vjust = -10, hjust = 0.1), strip.background = element_blank(), strip.text=element_text(hjust=-0.003)) +
  facet_wrap(~factor(Lake, levels = c("Frazer Lake", "Karluk Lake", "Red Lake")), ncol = 1)

#Create a high resolution TIFF file for publication
ggsave(filename = "Figure2_FrazerPaper.tiff", plot = Figure2, device = "tiff", height = 10, width = 6, units="in", 
       dpi=100)


# 3.5-Prepare M-ratio input file for microsatellite data from Frazer Lake (some modification of output is needed for use in M-ratio test)-----------------------
NewInput <- input_NucGoMsat_F
#Use 'separate' in dplyr package to split each genotype column into two columns of alleles for each locus 
for(i in 1:length(LociListNucGoMsat_F)) {
  NewInput <-  separate(NewInput, LociListNucGoMsat_F[i], into = c(paste(LociListNucGoMsat_F[i],"_1", sep = ""), paste(LociListNucGoMsat_F[i], "_2", sep = "")), sep = "/")  
}
#Vector of Frazer collection names
FrazerPops <- unique(NewInput$Pop)
#remove three microsatellite loci with complex repeats (not simple di- or tetra-nucleotide) 
NewInputS1 <- select(NewInput, -c(One_One100_1, One_One100_2, One_One105_1, One_One105_2, One_Ots107v4_1, One_Ots107v4_2))
#Number of loci used in M-ratio test
FrazerLociNo <- (length(NewInputS1)-2)/2
#Loop writes file for each population of allele counts for each locus 
for (i in 1:length(FrazerPops)){
  FrazerTemp <- filter(NewInputS1, Pop == FrazerPops[i])
  for (j in 1:FrazerLociNo){
    TempCount <- table(c(pull((FrazerTemp)[j*2+1]), pull((FrazerTemp)[j*2+2])), exclude = "000")
    write.table(t(TempCount), paste("FrazerNew", FrazerPops[i], ".txt", sep = ""),quote=FALSE,row.names=F,col.names=T,sep="\t",na="", append = TRUE)    
  }
}
#the next four lines are not needed
write.table(NewInput, "FrazerAlleles.txt", quote=FALSE,row.names=F,col.names=T,sep=" ",na="")
for (i in 1:length(LociListNucGoMsat_F)){
  print(table(c(pull((NewInput)[i*2+1]), pull((NewInput)[i*2+2])), exclude = "000")) #exludes "000"
}
#remove temp variables used for loop
rm("i", "j", "NewInput", "FrazerPops", "NewInputS1", "FrazerLociNo", "FrazerTemp", "TempCount")


# 3.6-Create Genepop input file from tibble using conversion function-----------------------------------------
#Genepop conversion function
Tib2Genepop <- function(x, y, z) {
  input2GP <- x
  #Optional: remove "One_' from prefix of each locus name (so can distinguish loci in output files from GenePop analysis)
  names(input2GP) <- gsub("One_","",names(input2GP),fixed=T) 
  #Creates blank dataframe with same number of columns as inputGP and no. columns - 1 rows.
  NewGP <-  as_tibble(matrix(nrow = length(input2GP) - 1, ncol = length(input2GP))) 
  #Assign the first column of NewGP the column names of inputGP column 2 through end.  
  NewGP[,1] <-  names(input2GP[,2:length(input2GP)]) 
  #Assign the column names of inputGP to NewGP
  names(NewGP) <-  names(input2GP) 
  #Assign the first row, column 1 of NewGP the string "Genepop Input'.
  NewGP[1,1] <- 'Genepop Input' 
  #Create LociCount vector with one element that is the number of columns of NewGP minus 2.
  LociCount <-  length(NewGP)-2 
  input2GP[,-c(1,2)] <- mutate_if(input2GP[,-c(1,2)], is.character, str_replace_all, pattern = "/", replacement = "")
  #Add a comma to each character in the column IndList in the dataframe inputGP.  The comma in front of the individual ID is needed for genepop. 
  input2GP$Ind <- paste(input2GP$Ind, ',', sep = '') 
  #Create a one row empty dataframe 
  Pop <-  as.data.frame(t(c(rep(NA, length(input2GP))))) 
  #Assign the column names of input2GP to Pop
  names(Pop) <-  names(input2GP) 
  #Assign the first row, column 1 of Pop the string 'Pop'
  Pop[1] <-  'Pop' 
  #Create a PopN vector with the names of each population from the the column input2GP$Pop in input2GP.
  PopN <-  unique(input2GP$Pop)
  #for loop appends the datafram NewGP with the Pop dataframe and each population, in turn, from the dataframe inputGP.
  for(i in 1:length(PopN)) {
    tempGP <- filter(input2GP, Pop == PopN[i])
    NewGP <-  bind_rows(NewGP, Pop, tempGP)  
  }
  #Remove the column titled PopList from NewGP
  NewGP <-  select(NewGP, -c(Pop)) 
  #Write NewGP to file.
  write.table(NewGP,y,quote=FALSE,row.names=F,col.names=F,sep=" ",na="") 
  write.table(PopN,z, quote=FALSE,row.names=F,col.names=F,sep=" ",na="")
}

#Here, calls the nuclear+mtDNA input file.  Can also use the nuclear loci file only.
Tib2Genepop(input_NucMtDNA_GoL, "Kodiak_NucmtDNA_GP.txt", "Kodiak_NucmtDNA_GP_PopN.txt")


# 4=SECOND LEVEL ANALYSIS---------------------------------------------------------------------------------------------

# 4.1-Pairwise G-test of genotypic frequency homogeniety among all collections-------------------------------------------

#Genepop
#G-test of differentiation - overall and pairwise for combined data (nuclear loci and mtDNA)
locinfile <-  "Kodiak_NucmtDNA_GP.txt"

#Overall
test_diff(locinfile, genic = F, pairs = F, outputFile = "Kodiak_Gtest_genepop.txt", dememorization = 10000, batches = 100, iterations = 5000)

#Pairwise
test_diff(locinfile, genic = F, pairs = T, outputFile = "Kodiak_GtestPW_genepop.txt", dememorization = 10000, batches = 100, iterations = 5000)

#Prepare PW G-test results for heat map (used "Kodiak_GtestPW_genepop_sub.txt" created by hand from "Kodiak_GtestPW_genepop.txt")
PW_Gtest <- read_table("Kodiak_GtestPW_genepop_sub.txt", skip = 5, col_names = F) %>% #parsing failure warning is due to the term "Highly Sign" in col 5 in some rows 
  select(X1, X3, X5) %>% 
  rename(Pop1 = X1, Pop2 = X3, GtestPval = X5) %>% 
  mutate(GtestPval = recode(GtestPval, "Highly" = "0.000000")) %>% 
  mutate(GtestPval = str_remove(GtestPval, "<"))
PW_Gtest_PopsCode <- unique(c(PW_Gtest$Pop1, PW_Gtest$Pop2))
PW_Gtest_PopNum <- as.character(c(1:17,19:33))
PW_Gtest_PopsCodeAbbr <- substr(PW_Gtest_PopsCode, 1,3)
PW_Gtest_PopsCodeAbbrNum <- paste(PW_Gtest_PopsCodeAbbr, PW_Gtest_PopNum, sep = "")
PW_Gtest_PopsCodeNum1 <- setNames(PW_Gtest_PopNum, PW_Gtest_PopsCode)
PW_Gtest_PopsCodeNum2 <- setNames(PW_Gtest_PopsCodeAbbrNum, PW_Gtest_PopsCode)
PW_Gtest$Pop1 <- recode(PW_Gtest$Pop1, !!!PW_Gtest_PopsCodeNum2)
PW_Gtest$Pop2 <- recode(PW_Gtest$Pop2, !!!PW_Gtest_PopsCodeNum2)
PW_Gtest$GtestPval <- as.numeric(PW_Gtest$GtestPval)
TempTib <- tibble(Pop1 = PW_Gtest_PopsCodeAbbrNum, Pop2 = PW_Gtest_PopsCodeAbbrNum, GtestPval = rep(1, 32))
#TempTib <- tibble(Pop1 = PW_Gtest_PopNum, Pop2 = PW_Gtest_PopsCodeAbbrNum, GtestPval = rep(1, 32))
PW_Gtest <- bind_rows(TempTib, PW_Gtest)
#PW_Gtest <- mutate(PW_Gtest, sig = ifelse(GtestPval < 0.05, "1", "0"))
PW_Gtest <- mutate(PW_Gtest, GtestSig = ifelse(GtestPval < 0.001, "3", ifelse(GtestPval < 0.01, "2", ifelse(GtestPval < 0.05, "1", 0))))
PopOrder <- unique(PW_Gtest$Pop1)
PopOrder2 <- unique(PW_Gtest$Pop2)
PW_Gtest$Pop1 <- factor(PW_Gtest$Pop1, levels = PopOrder)
PW_Gtest$Pop2 <- factor(PW_Gtest$Pop2, levels = PopOrder2)
PW_Gtest <- arrange(PW_Gtest, Pop2, Pop1) %>% 
  mutate(Pop1Pop2 = paste(Pop1,Pop2, sep=""))

SigTests <- count(PW_Gtest, PW_Pval <= 0.05)#not used in heatmap


# 4.2-Pairwise Fst among all collections---------------------------------------------------------------------

#Prepare data for use in adegenet, hierfsat packages
#Read all-loci data file "input_NucMtDNA_GoL"
input_NucMtDNA_GoL <- readRDS("input_NucMtDNA_GoL.rds")

#This code filters out all but the donor populations and Frazer Lake.  
#input_NucMtDNA_GoL <- filter(input_NucMtDNA_GoL, str_detect(Pop, "^F|KBL13|RTE31|RUL33"))
#This code filters out all but the Karluk and Red Lake donors and Frazer Lake. 
#input_NucMtDNA_GoL <- filter(input_NucMtDNA_GoL, str_detect(Pop, "^F|KBL13|RTE31"))
#This code filters out all but the microsatellite loci
#input_NucMtDNA_GoL <- select(input_NucMtDNA_GoL, Ind, Pop, One_Omy77v4:One_uSat60v4) %>% 
#  filter_at(vars(One_Omy77v4:One_uSat60v4), all_vars(. != "000/000"))
#This code filters out all but the SNP loci and removes indivdiduals with no score for all loci
#input_NucMtDNA_GoL <- select(input_NucMtDNA_GoL, Ind, Pop, One_agt_132:One_zP3b, HaploCode2) %>% 
#  filter_at(vars(One_agt_132:HaploCode2), all_vars(. != "000/000"))

#Vector of populations from nuclear data (for adegenet)
PopList_NucMtDNA_GoL <- input_NucMtDNA_GoL$Pop #List of populations for Adegent
#Vector of individuals from nuclear data (for adegenet)
IndList_NucMtDNA_GoL <- input_NucMtDNA_GoL$Ind #List of individuals from Adegent
#Vector of all loci
LociList_NucMtDNA_GoL <- names(select(input_NucMtDNA_GoL, -c(Pop, Ind)))
#Use 'select' in dplyr package to remove columns Ind, Pop.
input_NucMtDNA_GoL4AD <- select(input_NucMtDNA_GoL, -c(Ind, Pop)) 
#Convert to adegenet input file
AD_NucMtDNA_GoL <- df2genind(input_NucMtDNA_GoL4AD,sep = '/', ncode=6, ind.names=IndList_NucMtDNA_GoL, pop=PopList_NucMtDNA_GoL, NA.char="000/000", type="codom") 
#Convert to hierfstat input file
HF_NucMtDNA_GoL <- genind2hierfstat(AD_NucMtDNA_GoL) 
#Use Adegenet to convert genind object to genpop object for Poppr
Gen_NucMtDNA_GoL <- genind2genpop(AD_NucMtDNA_GoL)

#Use Hierfstat to compute pairwise Fst using all loci
WCpw <- pairwise.WCfst(HF_NucMtDNA_GoL) %>% 
  round(3) %>% 
  replace(is.na(.), 0)

#Convert pairwise FST matrix to column form for merging with pairwise G-test results
WCpw_melt = melt(WCpw, na.rm = TRUE) %>% 
  mutate(value = replace(value, value < 0, 0)) %>% 
  rename("Fst" = "value") %>% 
  rename("Pop1" = "Var1", "Pop2" = "Var2") %>% 
  mutate(Pop1Pop2 = paste(Pop1,Pop2, sep="")) #create variable column for joining with G-test results 


# 4.3-FIGURE 3: G-test results superimposed on heatmap of pairwise Fst for nuclear loci and mtDNA----------------------------- 

#Join data tables of Fst and G-test results
FstGtest <- left_join(WCpw_melt, PW_Gtest, by = c("Pop1Pop2"), keep = FALSE) %>% 
  rename("Pop1" = "Pop1.x", "Pop2" = "Pop2.x") %>% 
  select(Pop1, Pop2, Fst, GtestSig) %>% 
  #pivot_longer(-c(1,2), names_to = "testN", values_to = "val") %>% 
  mutate(Pop1n = as.integer(substr(Pop1,4,5)), Pop2n = as.integer(substr(Pop2,4,5)))

#this code is used to compute mean and range of Fst values within and among lakes
UniqFstGtest <- filter(FstGtest, Pop1n < Pop2n) %>% 
  mutate(Pop2 = str_replace_all(Pop2, c("RUL33" = "BEC33"))) 
DifLakesFstGtest <- filter(UniqFstGtest, substr(Pop1,1,1) != substr(Pop2,1,1))
SameLakesFstGtest <- filter(UniqFstGtest, substr(Pop1,1,1) == substr(Pop2,1,1))

#Plot Gtest results as points superimposed on heat map of pairwise Fst
Figure3 <- ggplot(mapping = aes(Pop1, Pop2)) +
  geom_tile(data = FstGtest %>% filter(Pop1n <= Pop2n), aes(fill = Fst), color = "gray20") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 0.0425, limit = c(0,0.085), space = "Lab", name="Fst") + 
  #geom_tile(data = test3 %>% filter(Pop1n < Pop2n), fill = "transparent", color = "gray") +
  geom_point(data = FstGtest %>% filter(Pop1n < Pop2n, GtestSig != 0), aes(color = GtestSig), size = 2) +
  #geom_point(data = test3 %>% filter(testN == "GtestPval", Pop1n < Pop2n), aes(color = val)) +
  #scale_color_gradient2(low = "red", high = "white", mid = "orange", midpoint = 0.01, limit = c(0,1), space = "Lab", name="GtestPval") + 
  #scale_color_manual(values = c("white", "lightblue", "dodgerblue", "darkblue"), breaks = c(0,1:3), labels = c(">=0.050", " <0.050", "<0.010", "<0.001"), name="GtestPval")  +
  scale_color_manual(values = c("lightblue", "dodgerblue", "darkblue"), breaks = c(1:3), labels = c("<0.050", "<0.010", "<0.001"), name="Gtest_P")  +
  theme_minimal() + 
  coord_fixed(clip = "off") +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.5, size = 5, hjust = 0.5), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, hjust = 0.5, color = "black"), 
        #axis.ticks = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        #axis.line = element_line(colour = "black", linewidth = 0.25), 
        #panel.border = element_rect(colour = "black", fill="transparent", linewidth=0.5), 
        #        axis.line = element_line(colour = "black"), 
        #        panel.border = element_rect(colour = "blue", fill=NA, linewidth=3), 
        legend.direction = "horizontal",
        legend.position = c(0.80, 0.85), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12),  
        legend.box.just = "center",
        legend.key.size = unit(2,"line"),
        legend.margin = margin(-0.2,0,0,0, unit="cm"),
        plot.margin = margin(1,1,1,1, 'cm'), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  #scale_x_discrete(position = 'top', labels = c(1:32), expand = c(0,0)) + 
  #scale_x_discrete(position = 'top') + 
  guides(fill = guide_colorbar(title.position="top", title.hjust = 0.5)) +
  guides(color = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(size = 4))) +
  scale_y_discrete(limits = rev(levels(factor(FstGtest$Pop2)))) +
  annotate("text", x = 1:32, y = 32:1, label = "-", size = 2) + 
  annotate("text", x = 1:32 + 0, y = c(32:1) + 1.5, label = unique(FstGtest$Pop1), size = 3.5, angle = 90) +
  #  annotate("text", x = inputM_melt$Pop1, y = inputM_melt$Pop2, label = inputM_melt$sig, size = 2) +  
  geom_segment(data = FstGtest, aes(x = 4, xend = 11, y = 32 + 0.5, yend = 25 + 0.5), 
               color = "black", linewidth = 0.5, linetype = 2) +
  #geom_segment(data = FstGtest, aes(x = 1 - 0.5, xend = 1 - 0.5, y = 22 + 0.5, yend = 32 + 0.5), 
               #color = "black", linewidth = 0.5, linetype = 1) +
  #geom_segment(data = FstGtest, aes(x = 1 - 0.5, xend = 10 + 0.5, y = 22 + 0.5, yend = 22 + 0.5), 
               #color = "black", linewidth = 0.5, linetype = 1) +
  #geom_rect(data = FstGtest, aes(xmin = 1 - 0.5, xmax = 10 + 0.5, ymin = 32 + 0.5, ymax = 22 + 0.5),
            #fill = "transparent", color = "black", linewidth = 0.25, linetype = 2) +
  geom_segment(data = FstGtest, aes(x = 13, xend = 26, y = 23 + 0.5, yend = 10 + 0.5), 
               color = "black", linewidth = 0.5, linetype = 2) +
  #geom_segment(data = FstGtest, aes(x = 11 - 0.5, xend = 11 - 0.5, y = 7 + 0.5, yend = 22 + 0.5), 
               #color = "black", linewidth = 0.5, linetype = 1) +
  #geom_segment(data = FstGtest, aes(x = 11 - 0.5, xend = 25 + 0.5, y = 7 + 0.5, yend = 7 + 0.5), 
               #color = "black", linewidth = 0.5, linetype = 1) +  
  #geom_rect(data = FstGtest, aes(xmin = 11 - 0.5, xmax = 25 + 0.5, ymin = 22 + 0.5, ymax = 7 + 0.5),
            #fill = "transparent", color = "black", linewidth = 0.25, linetype = 2) +
  geom_segment(data = FstGtest, aes(x = 28, xend = 32, y = 8 + 0.5, yend = 4 + 0.5), 
               color = "black", linewidth = 0.5, linetype = 2) +
  #geom_segment(data = FstGtest, aes(x = 26 - 0.5, xend = 26 - 0.5, y = 1 + 0.5, yend = 7 + 0.5), 
               #color = "black", linewidth = 0.5, linetype = 1) +
  #geom_segment(data = FstGtest, aes(x = 26 - 0.5, xend = 31 + 0.5, y = 1 + 0.5, yend = 1 + 0.5), 
               #color = "black", linewidth = 0.5, linetype = 1) +    
  #geom_rect(data = FstGtest, aes(xmin = 26 - 0.5, xmax = 31 + 0.5, ymin = 7 + 0.5, ymax = 1 + 0.5),
            #fill = "transparent", color = "black", linewidth = 0.25, linetype = 2) +
  annotate("text", x = 8.5 + 0.5, y = 30, label = "Frazer", size = 6) + 
  annotate("text", x = 22 + 0.5, y = 17, label = "Karluk", size = 6) +  
  annotate("text", x = 31 + 0.5, y = 7, label = "Red", size = 6)

#Create a high resolution TIFF file for publication
ggsave(filename = "Figure3_FrazerPaper.tiff", plot = Figure3, device = "tiff", height = 10, width = 10, units="in", 
       dpi=100)


# 4.4-FIGURE 4: Scatter plots of MDS of PW Fst with and without locus One_U1004_183--------------------------

#Collection number from PW Fst matrix 
CollectionNo <- substr(colnames(WCpw), 4,5)

#Ecotype from PW Fst matrix
LakeHabitatTime <- substr(colnames(WCpw), 1, 3)
LakeHabitatTime <- replace(LakeHabitatTime, LakeHabitatTime == "RUL", "UOL")
#Ecotype <- substr(colnames(WCpw), 1, 2)

#Multi-dimensional scaling of PW Fst estimates for all collections
MDSFst <- cmdscale(WCpw,eig=FALSE,k=3) %>% 
  as_tibble(.name_repair = "unique") %>% 
  rename(Dimension_1 = 1, Dimension_2 = 2, Dimension_3 = 3) %>%   
  mutate(CollectionNo = CollectionNo, Lake = substr(LakeHabitatTime, 1, 1), Habitat = substr(LakeHabitatTime, 2, 2), ) %>% 
  select(CollectionNo, Lake, Habitat, everything())

#Color for each ecotype in plot
#myColEcotype = c("dodgerblue4", 'dodgerblue', 'deepskyblue', 'green', 'springgreen3', 'palegreen4', 'darkgreen', "gold2", "gold4", "slateblue")

#myColEcotype = c("dodgerblue4", 'dodgerblue4', 'dodgerblue4', 'palegreen4', 'palegreen4', 'palegreen4', 'palegreen4', "gold2", "gold2", "slateblue")

#myColEcotype = c("dodgerblue4", 'dodgerblue4', 'dodgerblue4', 'palegreen4', 'palegreen4', 'palegreen4', "gold2", "gold2", "slateblue")

myColEcotype = c("dodgerblue4", 'palegreen4', "gold2", "slateblue") #Color for each lake instead of lake + ecotype

#Scatter plot of MDS of PW Fst estimates including locus One_U1004_183
Figure4a <- ggplot(MDSFst, aes(x = Dimension_1, y = Dimension_2, color = Lake, shape = Habitat)) + 
#Data visualization  
  geom_point(stroke = 1) +
#layout
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=10), axis.text = element_text(size=8), 
        legend.key = element_rect(linewidth = 0.5), legend.key.size = unit(0.1, "lines"), 
        legend.position = "bottom", legend.title = element_text(size = 7), 
        legend.text = element_text(size = 6), legend.spacing.y = unit(0.25, 'lines'),
        legend.background = element_rect(fill='transparent')) +
#overrides
  scale_color_manual(values = myColEcotype) +
  #scale_shape_manual(values = c(rep(19,3), rep(2,4), rep(5,3), 25)) +
  #scale_shape_manual(values = c(17, 15, 18, 2, 0, 5, 6, 2, 5, 5)) +
  #scale_shape_manual(values = c(0, 1, 2, 0, 1, 2, 6, 0, 2, 1)) +
  #scale_shape_manual(values = c(0, 1, 2, 0, 1, 2, 0, 2, 1)) +
  scale_shape_manual(values = c(0, 1, 2)) +
  scale_y_continuous(limits = c(-0.04,0.05), breaks = seq(-0.04,0.05,0.02), name = "Dimension 2") +
  scale_x_continuous(limits = c(-0.04,0.06), breaks = seq(-0.04,0.06,0.02), name = "Dimension 1") +
  #scale_y_continuous(limits = c(-0.02,0.02), breaks = seq(-0.02,0.02,0.01), name = "Dimension 2") +
  #scale_x_continuous(limits = c(-0.02,0.03), breaks = seq(-0.02,0.03,0.01), name = "Dimension 1") +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  #guides(color = guide_legend(override.aes = list(size = 1.0), nrow = 2, byrow = TRUE, title = "test")) +
  #guides(color = guide_legend(nrow = 2, byrow = TRUE, title = "test")) +
  #annotate('text', x = 0.019, y = 0.022, label = "Karluk, Late\nKBL,KOL,KTL", size = 2,  color = "black") +
  #annotate('text', x = 0.022, y = 0.005, label = "Karluk, Early\nKTE", size = 2,  color = "black") +
  #annotate('text', x = -0.017, y = 0.005, label = "Red, Late\nRBL", size = 2,  color = "black") +
  #annotate('text', x = -0.021, y = -0.028, label = "Red, Early\nRTE", size = 2,  color = "black") +
  #annotate('text', x = -0.013, y = -0.014, label = "Frazer\nFBM,FOM,FTM", size = 2,  color = "black") +
  annotate('text', x = 0.019, y = 0.022, label = "Karluk, Late", size = 2,  color = "black") +
  annotate('text', x = 0.023, y = 0.005, label = "Karluk, Early", size = 2,  color = "black") +
  annotate('text', x = -0.017, y = 0.004, label = "Red, Late", size = 2,  color = "black") +
  annotate('text', x = -0.021, y = -0.028, label = "Red, Early", size = 2,  color = "black") +
  annotate('text', x = -0.010, y = -0.011, label = "Frazer", size = 2,  color = "black") +
  annotate('text', x = 0.047, y = -0.039, label = "Ruth, Late", size = 2,  color = "black") +
  annotate('text', x = MDSFst$Dimension_1[13], y = MDSFst$Dimension_2[13] + 0.003, label = "13", size = 2,  color = "black") +
  annotate('text', x = MDSFst$Dimension_1[31], y = MDSFst$Dimension_2[31] + 0.003, label = "31", size = 2,  color = "black") +
  annotate('text', x = MDSFst$Dimension_1[32], y = MDSFst$Dimension_2[32] + 0.003, label = "33", size = 2,  color = "black") +
  annotate('text', x = -0.027, y = 0.05, label = "A", size = 4,  color = "black", parse = TRUE)

#Compute pairwise Fst without locus One_U1004_183
#Remove outlier locus One_U1004_183 from input
input_NucMtDNA_GoL <- select(input_NucMtDNA_GoL, -c(One_U1004_183))
#Repeat steps 4.2.1 and 4.2.2 above using input without One_U1004_183 
#Vector of populations from nuclear data (for adegenet)
PopList_NucMtDNA_GoL <- input_NucMtDNA_GoL$Pop #List of populations for Adegent
#Vector of individuals from nuclear data (for adegenet)
IndList_NucMtDNA_GoL <- input_NucMtDNA_GoL$Ind #List of individuals from Adegent
#Vector of all loci
LociList_NucMtDNA_GoL <- names(select(input_NucMtDNA_GoL, -c(Pop, Ind)))
#Use 'select' in dplyr package to remove columns Ind, Pop.
input_NucMtDNA_GoL4AD <- select(input_NucMtDNA_GoL, -c(Ind, Pop)) 
#Convert to adegenet input file
AD_NucMtDNA_GoL <- df2genind(input_NucMtDNA_GoL4AD,sep = '/', ncode=6, ind.names=IndList_NucMtDNA_GoL, pop=PopList_NucMtDNA_GoL, NA.char="000/000", type="codom") 
#Convert to hierfstat input file
HF_NucMtDNA_GoL <- genind2hierfstat(AD_NucMtDNA_GoL) 
#Use Adegenet to convert genind object to genpop object for Poppr
Gen_NucMtDNA_GoL <- genind2genpop(AD_NucMtDNA_GoL)

#Use Hierfstat to compute pairwise Fst without locus One_U1004_183
WCpw <- pairwise.WCfst(HF_NucMtDNA_GoL) %>% 
  round(3) %>% 
  replace(is.na(.), 0)

#Multi-dimensional scaling of PW Fst estimates without locus One_U1004_183 for all collections
MDSFstNo1004 <- cmdscale(WCpw,eig=FALSE,k=3) %>% 
  as_tibble(.name_repair = "unique") %>% 
  rename(Dimension_1 = 1, Dimension_2 = 2, Dimension_3 = 3) %>%   
  mutate(CollectionNo = CollectionNo, Lake = substr(LakeHabitatTime, 1, 1), Habitat = substr(LakeHabitatTime, 2, 2), ) %>% 
  select(CollectionNo, Lake, Habitat, everything())

#Scatter plot of MDS of PW Fst estimates without locus One_U1004_183
Figure4b <- ggplot(MDSFstNo1004, aes(x = Dimension_1, y = Dimension_2, color = Lake, shape = Habitat)) + 
#Data visualization  
  geom_point(stroke = 1) +
#layout
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=10), axis.text = element_text(size=8), 
        legend.key = element_rect(linewidth = 0.5), legend.key.size = unit(0.1, "lines"), 
        legend.position = "bottom", legend.title = element_text(size = 7), 
        legend.text = element_text(size = 6), legend.spacing.y = unit(0.25, 'lines'), 
        legend.background = element_rect(fill='transparent')) +
#overrides
  scale_color_manual(values = myColEcotype) +
  #scale_shape_manual(values = c(rep(19,3), rep(2,4), rep(5,3), 25)) +
  #scale_shape_manual(values = c(17, 15, 18, 2, 0, 5, 6, 2, 5, 5)) +
  #scale_shape_manual(values = c(0, 1, 2, 0, 1, 2, 6, 0, 2, 1)) +
  scale_shape_manual(values = c(0, 1, 2)) +
  #scale_y_reverse(limits = c(0.05,-0.04), breaks = seq(-0.04,0.05,0.01), name = "Dimension 2") +
  scale_y_continuous(limits = c(-0.04,0.05), breaks = seq(-0.04,0.05,0.02), name = "Dimension 2") +
  scale_x_continuous(limits = c(-0.04,0.06), breaks = seq(-0.04,0.06,0.02), name = "Dimension 1") +
  #scale_y_continuous(limits = c(-0.02,0.02), breaks = seq(-0.02,0.02,0.01), name = "Dimension 2") +
  #scale_x_continuous(limits = c(-0.02,0.03), breaks = seq(-0.02,0.03,0.01), name = "Dimension 1") +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  #guides(color = guide_legend(override.aes = list(size = 1.5), nrow = 1)) +
  #annotate('text', x = 0.016, y = -0.01, label = "Karluk, Late\nKBL,KOL,KTL", size = 2,  color = "black") +
  #annotate('text', x = 0.017, y = -0.003, label = "Karluk, Early\nKTE", size = 2,  color = "black") +
  #annotate('text', x = -0.007, y = 0.002, label = "Red, Late\nRBL", size = 2,  color = "black") +
  #annotate('text', x = -0.030, y = 0.011, label = "Red, Early\nRTE", size = 2,  color = "black") +
  #annotate('text', x = -0.014, y = -0.010, label = "Frazer\nFBM,FOM,FTM", size = 2,  color = "black") +
  annotate('text', x = 0.022, y = -0.01, label = "Karluk, Late", size = 2,  color = "black") +
  annotate('text', x = 0.023, y = -0.003, label = "Karluk, Early", size = 2,  color = "black") +
  annotate('text', x = -0.002, y = 0.002, label = "Red, Late", size = 2,  color = "black") +
  annotate('text', x = -0.030, y = 0.011, label = "Red, Early", size = 2,  color = "black") +
  annotate('text', x = -0.019, y = -0.002, label = "Frazer", size = 2,  color = "black") +
  annotate('text', x = 0.033, y = 0.046, label = "Ruth, Late", size = 2,  color = "black") +
  annotate('text', x = MDSFstNo1004$Dimension_1[13], y = MDSFstNo1004$Dimension_2[13] - 0.003, label = "13", size = 2,  color = "black") +
  annotate('text', x = MDSFstNo1004$Dimension_1[31], y = MDSFstNo1004$Dimension_2[31] + 0.003, label = "31", size = 2,  color = "black") +
  annotate('text', x = MDSFstNo1004$Dimension_1[32], y = MDSFstNo1004$Dimension_2[32] + 0.003, label = "32", size = 2,  color = "black") +
  #annotate('text', x = -0.027, y = 0.050, label = expression("F"[ST]), size = 4,  color = "black")
  annotate('text', x = -0.027, y = 0.050, label = "B", size = 4,  color = "black", parse = TRUE)

#Create a high resolution TIFF file for publication
Figure4 <- ggarrange(Figure4a, Figure4b, ncol = 1, nrow = 2, 
                          heights = c(8, 8, 1), common.legend = TRUE, legend = "bottom")
#Save figure to file
ggsave(filename = "Figure4_FrazerPaper.tiff", plot = Figure4, device = "tiff",
       height = 8, width = 4, units="in", dpi=200, bg="white")


# 4.5-Procrustes analysis and test comparing MDS results for Fst with and without locus One_U1004_183-----------------
proc <- procrustes(MDSFst[3:4], MDSFstNo1004[3:4], symmetric = FALSE)
#test comparing MDS results
prot <- protest(MDSFst[3:4], MDSFstNo1004[3:4], permutations = how(nperm = 1000))


# 4.6-AMOVA - First level hierarchical population structure using lake affiliation------------------------------

#use Arlequin outside R to perform AMOVA and test grouping hypotheses about hierarchical structure

#Used input file "input_NucMtDNA_GoL", converted to genepop file using code above, then used Arlequin to convert the resulting genepop file to Arlequin input
#Used python script "ArlequinPopMod.py to change population names (the Arlequin conversion gives the populations a generic name (e.g., "Population1"))
#G1: Two groups - Frazer beach with Karluk, Frazer tribs with Red.
#G2: Two groups - Frazer beach and trib with Red.
#G3: Two groups - Frazer beach and trib with Karluk.
#G4: Two groups - Frazer beach and trib as a third group

#Read csv file that contains AMOVA results for lake-level hierarchical structure
Kodiak_AMOVA <- read_csv("KodiakAMOVA.csv") %>% 
  #remove scenario G3 - Frazer collections grouped with Karluk - for manuscript
  filter(Strategy != "Scenario 3")

#Optional: create table for AMOVA results
TableAmova <- select(Kodiak_AMOVA, Strategy:Val) %>% 
  mutate_if(is.numeric, round, digits = 3) %>% 
  pivot_wider(names_from = Fstatistic, values_from = Val)

#Optional: Table of first level AMOVA results for four grouping scenarios for Frazer Lake collections relative to Red and Karluk lakes

write.csv(TableAmova, "Table2_FrazerPaper.csv", row.names = F)
write.table(TableAmova, "Table2_FrazerPaper.txt", quote=FALSE,row.names=F,col.names=T,sep=" ",na="")

# 4.7-FIGURE 5: Bar plot of first level AMOVA results-------------------------------------------------------------
#Did not include Group G3 in paper. Used G1, G2, G4 only.
StratNames <- c("Scenario 1\nRed+Fra(t)    Kar+Fra(b)", "Scenario 2\nRed+Fra    Kar", "Scenario 3\nFra   Red   Kar")
names(StratNames) <- c("Scenario 1", "Scenario 2", "Scenario 4")

Figure5 <- ggplot(Kodiak_AMOVA) +
  #Data visualization
  geom_bar(aes(x = Fstatistic, y = Val), stat = "identity", fill = "gray", color = "black") + 
  #Layout
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  #overrides
  scale_y_continuous(limits = c(0,0.035)) +
  scale_x_discrete(limits = c("Fst", "Fsc", "Fct")) +
  geom_text(aes(Fstatistic, Val, label=sprintf('%.3f',Val)), vjust=-0.25, size = 2.5) +
  geom_text(aes(Fstatistic, Val, label = "*"), vjust = -1.0, hjust = -3.5, size = 2) +
  #geom_errorbar(aes(x = Fstatistic, ymin=FctL95, ymax=FctU95), width = 0.15) +
  facet_wrap(~Strategy, ncol = 1, labeller = labeller(Strategy = StratNames))

#Print FIGURE 5 to directory
ggsave(filename = "Figure5_FrazerPaper.tiff", plot = Figure5, device = "tiff",
       height = 5.0, width = 3.5, units="in", dpi=200) #creates a high resolution TIFF file for publication


# 4.8-AMOVA - Second level hierarchical population structure using ecotype--------------------------------

input_NucMtDNA_GoL_Kar <- filter(input_NucMtDNA_GoL, str_sub(Pop,1,1) == "K")
input_NucMtDNA_GoL_Fra <- filter(input_NucMtDNA_GoL, str_sub(Pop,1,1) == "F")
input_NucMtDNA_GoL_Red <- filter(input_NucMtDNA_GoL, str_sub(Pop,1,1) == "R", Pop != "RUL1")

Tib2Genepop(input_NucMtDNA_GoL_Kar, "Kodiak_NucmtDNA_GP_Kar.txt", "Kodiak_NucmtDNA_GP_PopN_Kar.txt")
Tib2Genepop(input_NucMtDNA_GoL_Fra, "Kodiak_NucmtDNA_GP_Fra.txt", "Kodiak_NucmtDNA_GP_PopN_Fra.txt")
Tib2Genepop(input_NucMtDNA_GoL_Red, "Kodiak_NucmtDNA_GP_Red.txt", "Kodiak_NucmtDNA_GP_PopN_Red.txt")

input_Kar_NoU1004 <- select(input_NucMtDNA_GoL_Kar, -(One_U1004_183))
input_Fra_NoU1004 <- select(input_NucMtDNA_GoL_Fra, -(One_U1004_183))
input_Red_NoU1004 <- select(input_NucMtDNA_GoL_Red, -(One_U1004_183))

Tib2Genepop(input_Kar_NoU1004, "Kodiak_Kar_GP_NoU1004.txt", "Kodiak_Kar_GP_NoU1004_PopN.txt")
Tib2Genepop(input_Fra_NoU1004, "Kodiak_Fra_GP_NoU1004.txt", "Kodiak_Fra_GP_NoU1004_PopN.txt")
Tib2Genepop(input_Red_NoU1004, "Kodiak_Red_GP_NoU1004.txt", "Kodiak_Red_GP_NoU1004_PopN.txt")


# 4.9-FIGURE 6 Bar plot of second level AMOVA results----------------------------------------------------
# Read csv file that contains AMOVA results for ecotype-level hierarchical structure
Kodiak_AMOVA_InLakes <- read_csv("KodiakAMOVA_WithinLakes.csv") %>% 
  filter(U1004_YN == "w U1004_183")

Figure6 <- ggplot(Kodiak_AMOVA_InLakes) + 
#Data visualization
  geom_bar(aes(Fstatistic, Val), fill = "gray", color = "black", position = 'dodge',  stat = 'identity') + 
  #scale_fill_manual(values = c("black", "gray", "white")) +
#layout
  theme_bw() + 
  #ggtitle("(a)") +
#overides
  scale_y_continuous(name = "Val", breaks = seq(0,0.02,0.01)) +
  coord_cartesian(ylim = c(0,0.025)) +
  scale_x_discrete(limits = c("Fst", "Fsc", "Fct")) +
  geom_text(aes(Fstatistic, Val, label=sprintf('%.3f',Val)), vjust=-0.35, size = 2.5) +
  geom_text(aes(Fstatistic, Val, label = ifelse(Signif == "ns", "", Signif)), vjust = -1.0, hjust = -3.5, size = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=10), axis.text = element_text(size=8), 
        plot.title = element_text(size=12, vjust = 0, hjust = 0.05), legend.position = "none") +
  facet_wrap(~factor(Lake, levels = c("Frazer Lake ecotypes", "Karluk Lake ecotypes", "Red Lake ecotypes")), ncol = 1)

#creates a high resolution TIFF file for publication
ggsave(filename = "Figure6_FrazerPaper.tiff", plot = Figure6, device = "tiff",
       height = 5.0, width = 3.5, units="in", dpi=200) #creates a high resolution TIFF file for publication


# 4.10-DAPC using Adegenet----------------------------------------------------------------
#DAPC grouped by ecotypes, Ruth Lake removed
#Read all-loci file "input_NucMtDNA_GoL"
input_NucMtDNA_GoL <- readRDS("input_NucMtDNA_GoL.rds")
#Remove Ruth Lake collection
inputNucGoNoRuth <- filter(input_NucMtDNA_GoL, Pop != "RUL33")
#Use 'select' in dplyr package to remove columns Ind, Pop.
inputNucGoNoRuth4AD <- select(inputNucGoNoRuth, -c(Ind, Pop))
#Vector of populations from nuclear data (for adegenet)
PopListNucGoNoRuth <- inputNucGoNoRuth$Pop
#Vector of individuals from nuclear data (for adegenet)
IndListNucGoNoRuth <- inputNucGoNoRuth$Ind
#Convert to adegenet input file
ADnucGoNoRuth_DAPC <- df2genind(inputNucGoNoRuth4AD,sep = '/', ncode=6, ind.names=IndListNucGoNoRuth, pop=str_sub(PopListNucGoNoRuth,1,3), NA.char="000/000", type="codom") 
#Use ecotype code to group individuals for DAPC
ADnucGoNoRuth_DAPC_Ecotypes = ADnucGoNoRuth_DAPC$pop
#Alternative - use habitat code to group individuals for DAPC


#Optional code for doing correspondance analysis at the population level
#ForCA <- df2genind(inputNucGoNoRuth4AD,sep = '/', ncode=6, ind.names=IndListNucGoNoRuth, pop=PopListNucGoNoRuth, NA.char="000/000", type="codom") 
#ADnucGoNoRuth_CA <- genind2genpop(ForCA)
#CA_NoRuth <- dudi.coa(df = tab(ADnucGoNoRuth_CA), scannf=TRUE, nf=2)
#s.label(CA_NoRuth$li, sub="CA 1-2",csub=2)

#Cross validation for DAPC
xvalecotypes <- xvalDapc(tab(ADnucGoNoRuth_DAPC, NA.method = "zero"), ADnucGoNoRuth_DAPC_Ecotypes, training.set = 0.8)
xval2 <- xvalDapc(tab(ADnucGoNoRuth_DAPC, NA.method = "zero"), ADnucGoNoRuth_DAPC_Ecotypes, training.set = 0.8,
                  result = "groupMean", center = TRUE, scale = FALSE,
                  n.pca.max = 300, n.rep = 30, xval.plot = TRUE)

#DAPC - Based on results from cross validation - retained 100 PC's and 7 DF's
dapc1NoRuth<-dapc(ADnucGoNoRuth_DAPC,ADnucGoNoRuth_DAPC_Ecotypes, n.pca = 100, n.da = 7)
#Dataframe of LD coordinates 1 and 2 for each individual
dapc1NoRuthDF_12 = as_tibble(dapc1NoRuth$ind.coord[,1:2], rownames = "Ind") %>% 
  mutate(Ecotype = dapc1NoRuth$grp) %>% 
  mutate(Lake = substr(dapc1NoRuth$grp, 1, 1)) %>%
  mutate(Habitat = substr(dapc1NoRuth$grp, 2, 2)) %>% 
  mutate(LakeTime = paste(substr(dapc1NoRuth$grp, 1, 1), substr(dapc1NoRuth$grp, 3, 3), sep = "")) %>% 
  select(Ind, Ecotype, Lake, Habitat, LakeTime, everything())
#Dataframe of LD coordinates 1 and 2 for Frazer Lake individuals and Red and Karluk lake donors
dapc1NoRuthDF_12_Frazer <- filter(dapc1NoRuthDF_12, str_detect(Ind, "FTM|FBM|FOM|KBL4|RTE1"))
#Dataframe of LD coordinates 1 and 2 without Frazer Lake individuals
dapc1NoRuthDF_12_NoFrazer <- filter(dapc1NoRuthDF_12, !str_detect(Ind, "FTM|FBM|FOM"))
#Dataframe of LD coordinates 1 and 2 for each ecotype group
dapc1NoRuthDFgroup_12 = as_tibble(dapc1NoRuth$grp.coord[,1:2], rownames = "Ecotype") %>% 

  
#OPTIONAL - Dataframe of LD coordinates 2 and 3 for each individual
dapc1NoRuthDF_23 = as_tibble(dapc1NoRuth$ind.coord[,2:3], rownames = "Ind") %>% 
  mutate(Ecotype = dapc1NoRuth$grp) %>% 
  select(Ind, Ecotype, everything())
dapc1NoRuthDF_23_Frazer <- filter(dapc1NoRuthDF_23, substr(Ecotype, 1,1) == "F")

#Color for each ecotype
#myColNoRuthEcotype = c("dodgerblue4", 'dodgerblue', 'deepskyblue', 'green', 'springgreen3', 'palegreen4', 'darkgreen', "gold2", "gold4")
#myColNoRuthEcotype = c("dodgerblue4", "dodgerblue4", "dodgerblue4", 'palegreen4', 'palegreen4', 'palegreen4', 'palegreen4', "gold2", "gold2")
#Color for each lake (excluding Ruth)
myColNoRuthLake = c("dodgerblue4", 'palegreen4', "gold2") #Used instead of color for each lake + ecotype
#Line type for early and late time fish from Karluk and Red lakes
myColNoRuthLake_LakeTimeLines = c("dashed", "dotted", "dashed", "dotted")
#Label for plot legend
#GrpLabNoRuth_Ecotype = paste(c(1:9), '-', unique(dapc1NoRuthDF_12$Ecotype), sep = '')
#Group label for LD1 v LD2
GrpLabNoRuth_Ecotype <- unique(dapc1NoRuthDF_12$Ecotype)
#Group label for LD2 v LD3
GrpLabNoRuth_Ecotype23 <- unique(dapc1NoRuthDF_23$Ecotype)


#DAPC plot of LD1 and LD2 grouped by ecotype
dapc1_plot_12_ecotypes <- ggplot(dapc1NoRuthDF_12, aes(x = LD1, y = LD2, color = Ecotype, shape = Ecotype, size = Ecotype)) +
#Data visualization  
  geom_point() +
#layout
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_text(size = 10), axis.text.x = element_blank(), axis.title.x = element_text(size = 10), 
        legend.key = element_rect(linewidth = 0.5), legend.key.size = unit(0.01, "lines"), legend.position = c(0.90, 0.85), 
        legend.title = element_text(size = 6), legend.text = element_text(size = 5), legend.spacing.y = unit(0.15, 'lines'), 
        legend.background = element_rect(fill='transparent')) +
#overrides
  scale_color_manual(values = myColNoRuthEcotype, labels = GrpLabNoRuth_Ecotype) +
  #scale_shape_manual(values = c(rep(19,3), rep(2,4), rep(5,2)), labels = GrpLabNoRuth_Ecotype) +
  scale_shape_manual(values = c(17, 15, 18, 2, 0, 5, 6, 2, 5), labels = GrpLabNoRuth_Ecotype) +
  scale_size_manual(values = c(rep(0.75, 3), rep(0.25, 6)), labels = GrpLabNoRuth_Ecotype) +
  scale_y_continuous(limits = c(-4,4), name = "LD2") +
  scale_x_continuous(limits = c(-6,6), name = "LD1") +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  guides(color = guide_legend(override.aes = list(size = 1.25), byrow = TRUE)) +
  #stat_ellipse(data = dapc1NoRuthDF_12_Frazer, type = "norm", show.legend = F, size = 0.25, level = 0.95) +
  geom_mark_ellipse(data = dapc1NoRuthDF_12_Frazer, aes(color = Ecotype), size = 0.25) + #function from ggforce package, is not a statistical method, draws ellipse around all points
  annotate('label', x = dapc1NoRuth$grp.coord[,1], y = dapc1NoRuth$grp.coord[,2],  label = GrpLabNoRuth_Ecotype, size = 1.5, fill = myColNoRuthEcotype, color = "white") +
  annotate('label', x = -1.70, y = 0.09799645, label = "FBM", size = 1.5, fill = "dodgerblue4", color = "white") #added to adjust label to avoid overlap with FTM
  #annotate('text', x = dapc1NoRuth$grp.coord[,1], y = dapc1NoRuth$grp.coord[,2], label = c(1:9), size = 4,  color = "black")

#Alternative DAPC plot of LD1 and LD2 grouped by habitat
dapc1_plot_12_ecotypes <- ggplot(dapc1NoRuthDF_12, aes(x = LD1, y = LD2)) +
#Data visualization  
  geom_point(aes(colour = Lake, shape = Habitat)) +
#layout
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_text(size = 10), axis.text.x = element_blank(), axis.title.x = element_text(size = 10), 
        legend.key = element_rect(linewidth = 0.5), legend.key.size = unit(0.01, "lines"), legend.position = c(0.95, 0.80), 
        legend.title = element_text(size = 6), legend.text = element_text(size = 5), legend.spacing.y = unit(0.15, 'lines'), 
        legend.background = element_rect(fill='transparent')) +
#overrides
  #scale_color_manual(values = myColNoRuthLake, labels = GrpLabNoRuth_Ecotype) +
  scale_color_manual(values = myColNoRuthLake) +
  #scale_shape_manual(values = c(rep(19,3), rep(2,4), rep(5,2)), labels = GrpLabNoRuth_Ecotype) +
  #scale_shape_manual(values = c(17, 15, 18, 2, 0, 5, 6, 2, 5), labels = GrpLabNoRuth_Ecotype) +
  scale_shape_manual(values = c(0, 1, 2)) +
  #scale_size_manual(values = c(rep(0.75, 3), rep(0.25, 6)), labels = GrpLabNoRuth_Ecotype) +
  scale_size_manual(values = c(rep(0.75, 3))) +
  scale_linetype_manual(values = myColNoRuthLake_LakeTimeLines) +
  scale_y_continuous(limits = c(-4,4), name = "LD2") +
  scale_x_continuous(limits = c(-6,6), name = "LD1") +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  guides(color = guide_legend(override.aes = list(size = 2), byrow = TRUE)) +
  #stat_ellipse(data = dapc1NoRuthDF_12_NoFrazer, aes(group = LakeTime, color = Lake, linetype = LakeTime), type = "norm", show.legend = F, linewidth = 0.50, level = 0.95)
  geom_mark_ellipse(data = dapc1NoRuthDF_12, aes(group = LakeTime, color = Lake, filter = LakeTime != "FM", linetype = LakeTime, label = LakeTime), 
                    con.type = "none", label.fontsize = 8, label.fontface = c("plain"), label.buffer = unit(-3, "mm"), 
                    label.fill = F, show.legend = F, expand = 0.00005) #function from ggforce package, is not a statistical method, draws ellipse around all points
  #annotate('label', x = dapc1NoRuth$grp.coord[,1], y = dapc1NoRuth$grp.coord[,2],  label = GrpLabNoRuth_Ecotype, size = 1.5, fill = myColNoRuthEcotype, color = "white")
  #annotate('label', x = dapc1NoRuth$grp.coord[,1], y = dapc1NoRuth$grp.coord[,2],  label = substr(GrpLabNoRuth_Ecotype, 2,3), size = 3.0, fill = myColNoRuthEcotype, color = "white")
  #annotate('label', x = -1.70, y = 0.09799645, label = "FBM", size = 1.5, fill = "dodgerblue4", color = "white") #added to adjust label to avoid overlap with FTM
  #annotate('text', x = dapc1NoRuth$grp.coord[,1], y = dapc1NoRuth$grp.coord[,2], label = c(1:9), size = 4,  color = "black")


#Alternative DAPC plot of LD2 and LD3 grouped by ecotype
dapc1_plot_23_ecotypes <- ggplot(dapc1NoRuthDF_23, aes(x = LD2, y = LD3, color = Ecotype, shape = Ecotype, size = Ecotype)) + 
#Data visualization  
  geom_point() +
#layout
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_text(size = 10), axis.text.x = element_blank(), axis.title.x = element_text(size = 10), 
        legend.key = element_rect(linewidth = 0.5), legend.key.size = unit(0.01, "lines"), legend.position = c(0.90, 0.85), 
        legend.title = element_text(size = 6), legend.text = element_text(size = 5), legend.spacing.y = unit(0.15, 'lines'), 
        legend.background = element_rect(fill='transparent')) +
#overrides
  scale_color_manual(values = myColNoRuthEcotype, labels = GrpLabNoRuth_Ecotype23) +
  scale_shape_manual(values = c(rep(19,3), rep(2,4), rep(5,2)), labels = GrpLabNoRuth_Ecotype23) +
  scale_size_manual(values = c(rep(0.75, 3), rep(0.25, 6)), labels = GrpLabNoRuth_Ecotype23) +
  scale_y_continuous(limits = c(-4,4), name = "LD2") +
  scale_x_continuous(limits = c(-6,6), name = "LD1") +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  guides(color = guide_legend(override.aes = list(size = 1.25), byrow = TRUE)) +
  stat_ellipse(data = dapc1NoRuthDF_23_Frazer, type = "t", show.legend = F, size = 0.25, level = 0.95) +
  annotate('label', x = dapc1NoRuth$grp.coord[,2], y = dapc1NoRuth$grp.coord[,3],  label = GrpLabNoRuth_Ecotype23, size = 3, fill = myColNoRuthEcotype, color = "white")
  #annotate('text', x = dapc1NoRuth$grp.coord[,1], y = dapc1NoRuth$grp.coord[,2], label = c(1:9), size = 4,  color = "black")

#Histogram of DAPC Eigenvalues
EigenValues <- dapc1NoRuth$eig
dapcLD <- paste("LD", c(1:8), sep = "")
LD_EigenValues <- tibble(dapcLD, EigenValues)
EigenColor <- c(rep("gray20", 2), rep("gray", 6))

LD_eigenValues_plot <- ggplot(LD_EigenValues, aes(x = dapcLD, y = EigenValues)) +
#Data visualization
  geom_bar(stat = "identity", color = EigenColor, fill = EigenColor) +
#layout
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
#overrides
  scale_y_continuous(expand = c(0,0), limits = c(0, 1000)) +
  annotate('text', x = 5, y = 850, label = "DA eigenvalues", size = 2.5,  color = "black")

#DAPC plot of ecotypes with eigenvalue inset plot
dapc1_plot_12_ecotypes_eigenInset <- dapc1_plot_12_ecotypes + inset_element(LD_eigenValues_plot, left = 0.01, bottom = 0.01, right = 0.3, top = 0.2)


#DAPC plots of LD1 and LD2 loadings
#Dataframe of loadings for LD1 and LD2 by allele 
dapc1Load1 <- as_tibble(dapc1NoRuth$var.contr, rownames = "LocAllele") %>% 
  pivot_longer(cols = LD1:LD7, names_to = "Dimension", values_to = "Loading") %>% 
  filter(Dimension == "LD1" | Dimension == "LD2")

#DAPC LD loadings plots for dimensions 1 and 2
DAPC_LD12 <- ggplot(dapc1Load1) +
#Data visualization
  geom_bar(aes(x = LocAllele, y = Loading), stat = "identity", width = 0.25, color = "black") + 
#Layout
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 10), legend.position = "none") +
#overides
  geom_text_repel(data = filter(dapc1Load1, Loading > 0.10), 
                  aes(LocAllele, Loading, label = LocAllele), size = 2) +
  scale_y_continuous(limits = c(0,0.30), labels = scales::number_format(accuracy = 0.01)) +
  scale_x_discrete(limits = dapc1Load1$LocAllele, name = "Allele") +
  facet_wrap(~Dimension, ncol = 1)


# 4.11-FIGURE 7 : DAPC plot with LD1 and LD2 load plots----------------------------------------------------------------
#Combine DAPC plot and combined LD loading plots into single plot using ggarrange,
dapc1_Ecotypes_LD12_plot <- ggarrange(dapc1_plot_12_ecotypes_eigenInset,  DAPC_LD12, ncol = 1, nrow = 2,
                                      heights = c(4, 2))
#Print FIGURE 7 to directory
#dapc bar plot only
ggsave(filename = "Figure7_FrazerPaper.tiff", plot = dapc1_Ecotypes_LD12_plot, device = "tiff",
       height = 8, width = 6, units="in", dpi=125)


# 4.12-Estimate frequency of One_U1004_183 allele 1 for seasonal spawning ecotypes (early, mid, late) in each lake----------------------

U1004 <- select(input_NucMtDNA_GoL, Ind, Pop, One_U1004_183) %>% 
  separate(One_U1004_183, into = c('One_U1004_183_1', 'One_U1004_183_2'), sep = '/') %>% 
  gather('One_U1004_183_1', 'One_U1004_183_2', key = "LocusAllele", value = "Allele") %>% 
  filter(Allele != "000") %>%  
  mutate(PopGroup = substr(.$Pop, 1,3)) %>% 
  group_by(PopGroup, Pop)

U1004count <- count(U1004, Allele) %>% 
  mutate(freq = n/sum(n)) %>% 
  filter(Allele != "002") %>% 
  group_by(PopGroup) %>% 
  summarize(mean = mean(freq), lower = min(freq), upper = max(freq))

# 4.13-FIGURE S3 : plot of frequency of One_U1004_183 allele 1 for seasonal spawning ecotypes (early, mid, late) in each lake----------------------

U1004_plot <- ggplot(data = U1004count) + 
  geom_pointrange(mapping = aes(x=PopGroup, y=mean, ymin=lower, ymax=upper), size=0.5, stroke = 0.5, color="black", fill="white", shape=22) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_text(size = 7), 
        axis.text.y = element_text(size = 7), axis.title.y = element_text(size = 8)) +
  scale_y_continuous(name = "Freq", limits = c(0,1)) +
  scale_x_discrete(labels = c("FBM", "FOM", "FTM", "KBL", "KOL", "KTE", "KTL", "RBL", "RTE", "RUL")) + 
  annotate("text", x = c((1:10)-0.30), y = U1004count$mean, label = c(rep("Mid", 3), rep("Late", 2), "Early", rep("Late", 2), "Early", "Late"), size = 2) +
  annotate("text", x = 4, y = 1.0, label = "One_U1004_183 Allele 001", size = 3)

#creates a high resolution TIFF file for publication
ggsave('FigureS3_FrazerPaper.tiff', plot = U1004_plot, device = "tiff",
       height = 3.5, width = 5.0, units="in", dpi=200)


# 4.14-Test donor contribution to Frazer Lake collections using STRUCTURE---------------------------------------
#This analysis was done outside R using the STRUCTURE program
#This code creates a STRUCTURE input file and writes to project folder
PopFlagCon <- setNames(c(rep("0", 10), rep("1", 3)), c(1:10,13,31,33))
input_str <-   mutate(input_NucMtDNA_GoL, PopNo = str_sub(Pop, start = 4L, end = -1L)) %>% 
  select(Ind, Pop, PopNo, everything(), -c(HaploCode2)) %>% 
  filter(str_detect(Pop, "^F|KBL13|RTE31|RUL33")) %>% 
  pivot_longer(-c(1:3), names_to = "Loci", values_to = "Genotype") %>% 
  separate(Genotype, into = c("A1", "A2")) %>% 
  pivot_longer(c("A1", "A2"), names_to = "Allele", values_to = "Genotype") %>% 
  mutate(Genotype = str_replace_all(Genotype, c("000" = "-9"))) %>% 
  pivot_wider(names_from = "Loci", values_from = "Genotype") %>% 
  mutate(PopFlag = recode(PopNo, !!!PopFlagCon)) %>% 
  select(Ind, PopNo, PopFlag, Pop, everything())
#This code creates a string of the locus names that is appended as a header to the STRUCTURE file.
Loci <- t(names(input_str[-(1:5)]))
#Loci <- names(input_str)
#Append the STRUCTURE file to the file "Loci" and write to file.
write.table(Loci,"KodiakSock101_str.txt",quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ") 
write.table(input_str,"KodiakSock101_str.txt",quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", append = TRUE) 
#write.table(input_str,"KodiakSock101_str.txt",quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", Loci, append = TRUE) 

#The Structure output (10 interations = 10 files) was processed using the program Clumpp outside of R.
#The Clumpp output ("KodiakClumppRev") was read here. 
input_structure <- read.table("KodiakClumppRev") %>% 
  as_tibble() %>% 
  select(2,6:8) %>% 
  #rename columns - the order of donor populations is based on the Structure/clumpp output)
  rename("Ind" = "V2", "RTE" = "V6", "RUL" =  "V7", "KBL" = "V8") %>% 
  mutate("Pop" = substr(Ind, 1, 4)) %>% 
  mutate("Ecotype" = substr(Ind, 1, 3)) %>% 
  select("Ind", "Pop", "Ecotype", everything()) %>%
  filter(!(Ecotype %in% c("KBL", "RTE", "RUL"))) %>% 
  pivot_longer(-c(1,2,3), names_to = "Donor", values_to = "PostProb") 


ProportionSum <- group_by(input_structure, Ecotype, Donor) %>% 
  summarize(meanPostProb = mean(PostProb), varPostProb = var(PostProb)) %>% 
  pivot_wider(names_from = "Donor", values_from = "meanPostProb") %>% 
  mutate_if(is.double, round, 3)

#Vector of donor names to use as header for plot of individual admixture proportion data
DonorNames <- c("Red Lake Donor RTE31", "Ruth Lake Donor RUL33", "Karluk Lake Donor KBL13")
names(DonorNames) <- c("RTE", "RUL", "KBL")

# 4.15-FIGURE 8 : Boxplot of admixture proportions (Q) (from the three donor populations) for Frazer Lake samples grouped by ecotype----------------------------------------------------------------
#box-whisker plot of individual Q values - admixture proportions.
Figure8 <- ggplot(input_structure, aes(x = Ecotype, y = PostProb)) +
  #Data visualization  
  geom_boxplot(varwidth = TRUE, linewidth = 0.25) +
  #geom_jitter(width=0.05, alpha=0.5) +
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  stat_summary(fun = mean, geom="point", shape=15, size=2, fill="black") +
  #layout
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="gray")) +
  scale_y_continuous(limits = c(0,1.0)) +
  #overides
  ylab("Individual Admixture Proportion (Q)") +
  xlab("Frazer Lake Ecotype") +
  facet_grid(. ~ Donor, scales = "free_x", space = "free", labeller = labeller(Donor = DonorNames))

#Print FIGURE 8 to directory
#Boxplot
ggsave(filename = "Figure8_FrazerPaper.tiff", plot = Figure8, device = "tiff",
       height = 4, width = 6, units="in", dpi=200)

# 4.16-Alternative test of donor contribution to Frazer Lake collections Using DAPC and supplementary individuals-----------------------------------------------------------------
input_NucMtDNA_GoL <- readRDS("input_NucMtDNA_GoL.rds")
#Use 'select' in dplyr package to remove columns Ind, Pop.
input_NucMtDNA_GoL4AD <- select(input_NucMtDNA_GoL, -c(Ind, Pop))
#Vector of populations from nuclear data (for adegenet)
PopListinput_NucMtDNA_GoL <- input_NucMtDNA_GoL$Pop
#Vector of individuals from nuclear data (for adegenet)
IndListinput_NucMtDNA_GoL <- input_NucMtDNA_GoL$Ind
#Convert to adegenet input file
ADinput_NucMtDNA_GoL_DAPC <- df2genind(input_NucMtDNA_GoL4AD,sep = '/', ncode=6, ind.names=IndListinput_NucMtDNA_GoL, pop=str_sub(PopListinput_NucMtDNA_GoL,1,3), NA.char="000/000", type="codom") 
#Vector of collections representing the donor populations to use as training data
Donors <- subset(input_NucMtDNA_GoL$Ind, substr(input_NucMtDNA_GoL$Ind, 1, 4) %in% c("KBL4", "RTE1", "RUL1"))
#Remove some donors from Donor vector to use as test samples in self assign test.  Code works but did not use.
#DonorsRemoved <- sample(Donors, 60, replace = FALSE)
#DonorsRemain <- setdiff(Donors, DonorsRemoved)
#DAPC input for training data that include only Donors [Donors].  Could use [DonorsRemain] inplace of [Donors] if including a self assign test.  
Donors4DAPC <- ADinput_NucMtDNA_GoL_DAPC[Donors]
#Vector of ecotypes to use as group prior in DAPC
Donors4DAPC_Ecotypes <- substr(Donors4DAPC$pop, 1, 3)
#Vector of collections representing Frazer Lake as supplementary individuals
Frazer <- subset(input_NucMtDNA_GoL$Ind, substr(input_NucMtDNA_GoL$Ind, 1, 3) %in% c("FBM", "FOM", "FTM"))
#If using [DonorsRemain] for Doners4DAPC use this code for Frazer that adds DonorsRemoved for self assign test 
#Frazer <- subset(input_NucMtDNA_GoL$Ind, substr(input_NucMtDNA_GoL$Ind, 1, 3) %in% c("FBM", "FOM", "FTM")) %>% 
#c(DonorsRemoved)
#DAPC input for supplementary individuals.
Frazer4DAPC <- ADinput_NucMtDNA_GoL_DAPC[Frazer]
#Use cross validation to determine the optimum number of PCs to retain for DAPC
xval1 <- xvalDapc(tab(Donors4DAPC, NA.method = "zero"), Donors4DAPC_Ecotypes, training.set = 0.8)
#DAPC on training data [Donors] using 20 PCs based on cross validation results and retaining 2 axes.
dapcDonors <- dapc(Donors4DAPC, Donors4DAPC_Ecotypes, n.pca = 20, n.da = 2)
#Show scatter plot of DAPC results for training data samples [Donors]
scatter(dapcDonors)
#Perform DAPC on the supplementary samples [Frazer] using the DAPC results from the training samples [Donors]
predict_Frazer <- predict.dapc(dapcDonors, newdata = Frazer4DAPC)
#Tibble of posterior probability estimates for supplementary samples [Frazer] for plotting.
FrazerPosterior <- as_tibble(round(predict_Frazer$posterior,3), rownames = "Ind") %>%
  mutate(Pop = substr(Ind, 1, 4), Ecotype = substr(Ind, 1, 3)) %>% 
  select(Ind, Pop, Ecotype, everything()) %>% 
  pivot_longer(-c(1,2,3), names_to = "Donor", values_to = "PostProb")
#Summary of posterior probability estimates for supplementary samples [Frazer].
ProportionSum <- group_by(FrazerPosterior, Pop, Donor) %>% 
  summarize(meanPostProb = mean(PostProb)) %>% 
  pivot_wider(names_from = "Donor", values_from = "meanPostProb") %>% 
  mutate_if(is.double, round, 3)

#optional code that tests self assignment of Donor populations using leave two out.
FrazerDonorPosterior <- as_tibble(round(predict_Frazer$posterior,3), rownames = "Ind")
d = 1
for(i in 1:143) {
Donors <- subset(input_NucMtDNA_GoL$Ind, substr(input_NucMtDNA_GoL$Ind, 1, 4) %in% c("KBL4", "RTE1", "RUL1"))
DonorsRemoved <- Donors[d:(d+1)]
DonorsRemain <- setdiff(Donors, DonorsRemoved)
DonorsRemain4DAPC <- ADinput_NucMtDNA_GoL_DAPC[DonorsRemain]
DonorsRemain4DAPC_Ecotypes <- substr(DonorsRemain4DAPC$pop, 1, 3)
DonorsRemoved4DAPC <- ADinput_NucMtDNA_GoL_DAPC[DonorsRemoved]
dapcDonorsRemain <- dapc(DonorsRemain4DAPC, DonorsRemain4DAPC_Ecotypes, n.pca = 20, n.da = 2)
predict_DonorsRemoved <- predict.dapc(dapcDonorsRemain, newdata = DonorsRemoved4DAPC)
tempPosterior <- as_tibble(round(predict_DonorsRemoved$posterior,3), rownames = "Ind")
FrazerDonorPosterior <- bind_rows(FrazerDonorPosterior, tempPosterior)
d = d+2
}
#Tibble of posterior probability estimates for supplementary samples [Frazer and Donors] for plotting (donor estimates computed using leave two out).
FrazerDonorPosterior <- FrazerDonorPosterior %>% 
  mutate(Pop = substr(Ind, 1, 4), Ecotype = substr(Ind, 1, 3)) %>% 
  select(Ind, Pop, Ecotype, everything()) %>% 
  pivot_longer(-c(1,2,3), names_to = "Donor", values_to = "PostProb") 
#Summary of posterior probability estimates for supplementary samples [Frazer and Donors].
ProportionSumFrazerDonor <- group_by(FrazerDonorPosterior, Pop, Donor) %>% 
  summarize(meanPostProb = mean(PostProb)) %>% 
  pivot_wider(names_from = "Donor", values_from = "meanPostProb") %>% 
  mutate_if(is.double, round, 3)

#Individual Donor population samples with posterior probability estimates greater than 0.50 for the other population.
filter(FrazerDonorPosterior, Ecotype == "KBL", Donor == "RTE", PostProb > 0.50)
filter(FrazerDonorPosterior, Ecotype == "RTE", Donor == "KBL", PostProb > 0.50)

#Individual Frazer samples with posterior probability estimates greater than 0.50 for the Karluk Lake donor.
filter(FrazerPosterior, Ecotype %in% c("FBM", "FOM", "FTM"), Donor == "KBL", PostProb > 0.50)

# 4.17- Optional FIGURE 8 : Boxplot of posterior probability of membership in the three donor collections [Donor] for the supplementary samples [Frazer] grouped by ecotype------------------
DonorNames <- c("Karluk Lake Donor KBL13", "Red Lake Donor RTE31", "Ruth Lake Donor RUL33")
names(DonorNames) <- c("KBL", "RTE", "RUL")
DonorMembership_plot <- ggplot(FrazerPosterior, aes(x = Ecotype, y = PostProb)) +
  #Data visualization  
  geom_boxplot(varwidth = TRUE) +
  #geom_jitter(width=0.05, alpha=0.5) +
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  stat_summary(fun = mean, geom="point", shape=15, size=4, fill="black") +
  #layout
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="gray")) +
  #overides
  ylab("Posterior Membership Probability") +
  xlab("Frazer Lake Ecotype") +
  facet_grid(. ~ Donor, scales = "free_x", space = "free", labeller = labeller(Donor = DonorNames))

#Print FIGURE 8 to directory
#Boxplot
ggsave(filename = "Figure8_FrazerPaper.tiff", plot = DonorMembership_plot, device = "tiff",
       height = 3, width = 6, units="in", dpi=600)


# 4.18-Compare Hs and Ar in Frazer beach and tributary collections versus putative donor collections KBL4(Map no 13, beach), RTE1(Map no 31, tributary)-------------
#Hierfstat and Wilcoxon paired-sample test in R stats package

FrazerBeachPops <- unique(PopListNucGo)[1:4]
FrazerTribPops <- unique(PopListNucGo)[5:10]
FrazerBeachPopsTibble <- as_tibble(FrazerBeachPops)
FrazerTribPopsTibble <- as_tibble(FrazerTribPops)
mSatHs <- c()
mSatAr <- c()
SNPsHs <- c()
SNPsAr <- c()

B_mSatsGo <- basic.stats(HFnucGoMsat) 
AR_mSatsGo <- allelic.richness(HFnucGoMsat) 
B_SNPsGo <- basic.stats(HFnucGoSNP) 
AR_SNPsGo <- allelic.richness(HFnucGoSNP)
B_mtDNAGo <- basic.stats(HF_mtDNAGo, diploid = F)

#Beach - microsatellites - Hs
for(i in 1:length(FrazerBeachPops)) {
  tempHsWil <- wilcox.test(B_mSatsGo$Hs[,"KBL4"], B_mSatsGo$Hs[,FrazerBeachPops[i]], exact = F, paired = T, alternative = "g")
  mSatHs <-  c(mSatHs, tempHsWil$p.value)  
}
FrazerBeachPopsTibble <- mutate(FrazerBeachPopsTibble, mSatHs = mSatHs)
mSatHs <- c()

#Beach - microsatellites - Ar
for(i in 6:(length(FrazerBeachPops)+5)) {
  tempArWil <- wilcox.test(AR_mSatsGo$Ar[,13], AR_mSatsGo$Ar[,i], exact = F, paired = T, alternative = "g")
  mSatAr <-  c(mSatAr, tempArWil$p.value)  
}
FrazerBeachPopsTibble <- mutate(FrazerBeachPopsTibble, mSatAr = mSatAr)
mSatAr <- c()

#Beach - SNPs - Hs
for(i in 1:length(FrazerBeachPops)) {
  tempHsWil <- wilcox.test(B_SNPsGo$Hs[,"KBL4"], B_SNPsGo$Hs[,FrazerBeachPops[i]], exact = F, paired = T, alternative = "g")
  SNPsHs <-  c(SNPsHs, tempHsWil$p.value)  
}
FrazerBeachPopsTibble <- mutate(FrazerBeachPopsTibble, SNPsHs = SNPsHs)
SNPsHs <- c()

#Beach - SNPs - Ar
for(i in 6:(length(FrazerBeachPops)+5)) {
  tempArWil <- wilcox.test(AR_SNPsGo$Ar[,13], AR_SNPsGo$Ar[,i], exact = F, paired = T, alternative = "g")
  SNPsAr <-  c(SNPsAr, tempArWil$p.value)  
}
FrazerBeachPopsTibble <- mutate(FrazerBeachPopsTibble, SNPsAr = SNPsAr)
SNPsAr <- c()

#Trib - microsatellites - Hs
for(i in 1:length(FrazerTribPops)) {
  tempHsWil <- wilcox.test(B_mSatsGo$Hs[,"RTE1"], B_mSatsGo$Hs[,FrazerTribPops[i]], exact = F, paired = T, alternative = "g")
  mSatHs <-  c(mSatHs, tempHsWil$p.value)  
}
FrazerTribPopsTibble <- mutate(FrazerTribPopsTibble, mSatHs = mSatHs)
mSatHs <- c()

#Trib - microsatellites - Ar
for(i in 1:(length(FrazerTribPops))) {
  tempArWil <- wilcox.test(AR_mSatsGo$Ar[,31], AR_mSatsGo$Ar[,i], exact = F, paired = T, alternative = "g")
  mSatAr <-  c(mSatAr, tempArWil$p.value)  
}
FrazerTribPopsTibble <- mutate(FrazerTribPopsTibble, mSatAr = mSatAr)
mSatAr <- c()

#Trib - SNPs - Hs
for(i in 1:length(FrazerTribPops)) {
  tempHsWil <- wilcox.test(B_SNPsGo$Hs[,"RTE1"], B_SNPsGo$Hs[,FrazerTribPops[i]], exact = F, paired = T, alternative = "g")
  SNPsHs <-  c(SNPsHs, tempHsWil$p.value)  
}
FrazerTribPopsTibble <- mutate(FrazerTribPopsTibble, SNPsHs = SNPsHs)
SNPsAr <- c()

#Trib - SNPs - Ar
for(i in 1:(length(FrazerTribPops))) {
  tempArWil <- wilcox.test(AR_SNPsGo$Ar[,31], AR_SNPsGo$Ar[,i], exact = F, paired = T, alternative = "g")
  SNPsAr <-  c(SNPsAr, tempArWil$p.value)  
}
FrazerTribPopsTibble <- mutate(FrazerTribPopsTibble, SNPsAr = SNPsAr)
SNPsAr <- c()

# 4.19-TABLE X : Print each to output to screen and type results into Table 2----------------------------------------
write.csv(FrazerBeachPopsTibble, "TableX_FrazerPaper.csv")




# 5-UNUSED CODE ------------------------------------------------------------------------------
#works but not needed for paper
#DAPC grouped by Lakes, Ruth Lake included-----------------------------------------------------------------------------------------------
#Prepare for use in adegenet, hierfsat packages
#Vector of populations from nuclear data (for adegenet)
PopList_NucMtDNA_GoL <- input_NucMtDNA_GoL$Pop #List of populations for Adegent
PopList_NucMtDNA_GoL_noRuth <- filter(input_NucMtDNA_GoL$Pop == "RUL33") #List of populations for Adegent
#Vector of individuals from nuclear data (for adegenet)
IndList_NucMtDNA_GoL <- input_NucMtDNA_GoL$Ind #List of individuals from Adegent
#Vector of nuclear loci
LociList_NucMtDNA_GoL <- names(select(input_NucMtDNA_GoL, -c(Pop, Ind)))
#Use 'select' in dplyr package to remove columns Ind, Pop.
input_NucMtDNA_GoL4AD <- select(input_NucMtDNA_GoL, -c(Ind, Pop)) 
#Convert to adegenet input file
AD_NucMtDNA_GoL <- df2genind(input_NucMtDNA_GoL4AD,sep = '/', ncode=6, ind.names=IndList_NucMtDNA_GoL, pop=PopList_NucMtDNA_GoL, NA.char="000/000", type="codom") 
#Convert to hierfstat input file
HF_NucMtDNA_GoL <- genind2hierfstat(AD_NucMtDNA_GoL)

AD_NucMtDNA_GoL_Lakes = str_sub(AD_NucMtDNA_GoL$pop,1,3) %>% 
  str_replace_all(c("FBM|FOM|FTM" = "Frazer", "RBL|RTE" = "Red", "KBL|KTL|KTE|KOL" = "Karluk", "RUL" = "Ruth"))

dapc1 <- dapc(AD_NucMtDNA_GoL, AD_NucMtDNA_GoL_Lakes)

dapc1DF_12 <- as_tibble(dapc1$ind.coord[,1:2], rownames = "Ind") %>% 
  mutate(Lake = dapc1$grp) %>% 
  select(Ind, Lake, everything())

myColLakes = c("blue",  "green4",  'gray50', 'darksalmon')

GrpLabLakes = unique(dapc1DF_12$Lake)

# DAPC plot of PC1 and PC2
dapc1_plot_12 <- ggplot(dapc1DF_12, aes(x = LD1, y = LD2, color = Lake, shape = Lake)) + 
  geom_point(size = 0.5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(),
        axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), legend.key = element_rect(size = 0.5),
        legend.key.size = unit(0.01, "lines"), legend.position = c(0.90, 0.85), legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5), legend.spacing.y = unit(0.15, 'lines'), legend.background = element_rect(fill='transparent')) +
  scale_color_manual(values = myColLakes, labels = GrpLabLakes) +
  scale_shape_manual(values = c(21, 24, 23, 22), labels = GrpLabLakes) +
  scale_y_continuous(limits = c(-8,6)) +
  scale_x_continuous(limits = c(-6,6)) +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  guides(color = guide_legend(override.aes = list(size = 1.25), byrow = TRUE)) +
  #annotate('label', x = dapc1$grp.coord[,1], y = dapc1$grp.coord[,2], label = c(1:9), cex = 2, fill = myCol, color = "black") +
  stat_ellipse(type = "t", show.legend = F, size = 0.75) +
  #annotate('label', x = dapc1$grp.coord[,1], y = dapc1$grp.coord[,2], label = GrpLabLakes, cex = 1, fill = myColLakes, color = "white")
  annotate('text', x = -5.5, y = 5.5, label = "A", size = 4,  color = "black")

# DAPC plot of PC2 and PC3
dapc1DF_23 = as.data.frame(dapc1$ind.coord[,c(2,3)])
dapc1DF_23$Group = dapc1$grp
dapc1GroupTib_23 <- as_tibble(dapc1$grp.coord[,c(2,3)], rownames = "LakeGroup")
myCol = c('firebrick1', 'Firebrick2', 'firebrick3', 'green', 'green2', 'green4', 'darkgreen', 'dodgerblue', 'deepskyblue', 'darksalmon')
GrpLab = paste(c(1:10), '-', unique(dapc1DF$Group), sep = '')

dapc1_plot_23 <- ggplot(dapc1DF_23, aes(x = LD2, y = LD3, color = Group, shape = Group)) + 
  geom_point() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(),
        axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), legend.key = element_rect(size = 0.5),
        legend.key.size = unit(0.01, "lines"), legend.position = c(0.95, 0.70), legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5)) +
  scale_color_manual(values = myCol, labels = GrpLab) +
  scale_shape_manual(values = c(rep(1,3), rep(2,4), rep(5,2), 2), labels = GrpLab) +
  scale_y_continuous(limits = c(-8,7)) +
  scale_x_continuous(limits = c(-7,7)) +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  guides(color = guide_legend(override.aes = list(size = 1.75))) +
  #annotate('label', x = dapc1$grp.coord[,1], y = dapc1$grp.coord[,2], label = c(1:9), cex = 2, fill = myCol, color = "black") +
  #stat_ellipse(type = "t", show.legend = F, size = 0.75) +
  annotate('label', x = dapc1$grp.coord[,2], y = dapc1$grp.coord[,3], label = c(1:10), cex = 4, fill = myCol, color = "white")


#plot outlier loci results---------------------------------------------------------
FraSNPplot <- ggplot(OLinputF) + 
  geom_point(data = filter(OLinputF, FstPval > 0.01 & LocusName != c("MHC2_190", "MHC2_251")), aes(ObsH_BP, ObsFst), shape = 1) + 
  geom_point(data = filter(OLinputF, LocusName == c("MHC2_190", "MHC2_251")), aes(ObsH_BP, ObsFst), shape = 0, size = 3) +
  geom_point(data = filter(OLinputF, FstPval < 0.01), aes(ObsH_BP, ObsFst), shape = 2, size = 3) +
  theme_bw() +
  ggtitle("Frazer Lake") +
  scale_y_continuous(name = "", breaks = seq(0,0.4,0.1)) + 
  scale_x_continuous(name = "Ho/(1-Fst)", breaks = seq(0,0.6,0.1)) +
  geom_line(aes(Het_BP, `0.01`), linetype = "dashed", data = OLSimInputF) + 
  geom_line(aes(Het_BP, `0.99`), linetype = "dashed", data = OLSimInputF) + 
  coord_cartesian(xlim = c(0, 0.6), ylim = c(-0.05,0.40)) +
  geom_text(data = filter(OLinputF, LocusName == c("MHC2_190")), 
            aes(ObsH_BP, ObsFst, label = LocusName), nudge_x = -0.055, nudge_y = -0.015, size = 3) + 
  geom_text(data = filter(OLinputF, LocusName == c("MHC2_251")), 
            aes(ObsH_BP, ObsFst, label = LocusName), nudge_x = 0.075, nudge_y = -0.015, size = 3) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin=unit(c(0,0,0,0.5),"cm"), 
        axis.title = element_text(size=14), axis.text = element_text(size=12), 
        plot.title = element_text(size=14, vjust = -10, hjust = 0.1))

RedSNPplot <- ggplot(OLinputR) + 
  geom_point(data = filter(OLinputR, FstPval > 0.01 & LocusName != c("MHC2_190", "MHC2_251")), aes(ObsH_BP, ObsFst), shape = 1) + 
  geom_point(data = filter(OLinputR, LocusName == c("MHC2_190", "MHC2_251")), aes(ObsH_BP, ObsFst), shape = 0, size = 3) +
  geom_point(data = filter(OLinputR, FstPval < 0.01 & ObsH_BP > 0), aes(ObsH_BP, ObsFst), shape = 2, size = 3) +
  theme_bw() +
  ggtitle("Red Lake") +
  scale_y_continuous(name = "Fst", breaks = seq(0,0.4,0.1)) + 
  scale_x_continuous(name = "Ho/(1-Fst)", breaks = seq(0,0.6,0.1)) +
  geom_line(aes(Het_BP, `0.01`), linetype = "dashed", data = OLSimInputR) + 
  geom_line(aes(Het_BP, `0.99`), linetype = "dashed", data = OLSimInputR) + 
  coord_cartesian(xlim = c(0, 0.6), ylim = c(-0.05,0.40)) +
  geom_text(data = filter(OLinputR, LocusName == c("MHC2_190")), 
            aes(ObsH_BP, ObsFst, label = LocusName), nudge_x = -0.075, nudge_y = 0.0, size = 3) + 
  geom_text(data = filter(OLinputR, LocusName == c("MHC2_251")), 
            aes(ObsH_BP, ObsFst, label = LocusName), nudge_x = 0.075, nudge_y = 0.0, size = 3) + 
  geom_text(data = filter(OLinputR, FstPval < 0.01 & ObsH_BP > 0), 
            aes(ObsH_BP, ObsFst, label = LocusName), nudge_x = -0.075, nudge_y = 0.02, size = 3) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin=unit(c(0,0,0,0.5),"cm"), 
        axis.title = element_text(size=14), axis.title.x = element_blank(), axis.text = element_text(size=12), 
        plot.title = element_text(size=14, vjust = -10, hjust = 0.1))

KarSNPplot <- ggplot(OLinputK) + 
  geom_point(data = filter(OLinputK, FstPval > 0.01 & LocusName != c("MHC2_190", "MHC2_251")), aes(ObsH_BP, ObsFst), shape = 1) + 
  geom_point(data = filter(OLinputK, LocusName == c("MHC2_190", "MHC2_251")), aes(ObsH_BP, ObsFst), shape = 0, size = 3) +
  geom_point(data = filter(OLinputK, FstPval < 0.01), aes(ObsH_BP, ObsFst), shape = 2, size = 3) +
  theme_bw() +
  ggtitle("Karluk Lake") +
  scale_y_continuous(name = "", breaks = seq(0,0.4,0.1)) + 
  scale_x_continuous(name = "Ho/(1-Fst)", breaks = seq(0,0.6,0.1)) +
  geom_line(aes(Het_BP, `0.01`), linetype = "dashed", data = OLSimInputK) + 
  geom_line(aes(Het_BP, `0.99`), linetype = "dashed", data = OLSimInputK) + 
  coord_cartesian(xlim = c(0, 0.6), ylim = c(-0.05,0.40)) +
  geom_text(data = filter(OLinputK, LocusName == c("MHC2_190")), 
            aes(ObsH_BP, ObsFst, label = LocusName), nudge_x = 0.075, nudge_y = -0.015, size = 3) + 
  geom_text(data = filter(OLinputK, LocusName == c("MHC2_251")), 
            aes(ObsH_BP, ObsFst, label = LocusName), nudge_x = -0.075, nudge_y = 0.01, size = 3) + 
  geom_text(data = filter(OLinputK, FstPval < 0.01), 
            aes(ObsH_BP, ObsFst, label = LocusName), nudge_x = -0.075, nudge_y = 0.02, size = 3) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin=unit(c(0,0,0,0.5),"cm"), 
        axis.title = element_text(size=14), axis.title.x = element_blank(), 
        axis.text = element_text(size=12), 
        plot.title = element_text(size=14, vjust = -10, hjust = 0.1))

#Unused code for plotting outlier loci
#SNPoutliers for Red and Karluk lakes
OLinputRK <- read_tsv("fdist2_ObsOut_RK.txt") #read input as a tibble.
OLinputRK <- select(OLinputRK, c(1:5))
names(OLinputRK) <- c("LocusNo", "ObsH_BP", "ObsFst", "FstPval", "One_FstQuant")
OLinputRK <- bind_cols(OLinputRK, as.tibble(gsub("One_","",LociListNucGoSNP_RK,fixed=T)))
OLinputRK <- rename(OLinputRK, LocusName = value)
OLinputRK <- select(OLinputRK, LocusNo, LocusName, everything())
OLSimInputRK <- read_table("FstvHoNull_SNP_RK.txt")

ggplot(OLinputRK) + 
  geom_point(data = filter(OLinputRK, FstPval > 0.01 & LocusName != c("MHC2_190", "MHC2_251")), aes(ObsH_BP, ObsFst)) + 
  geom_point(data = filter(OLinputRK, LocusName == c("MHC2_190", "MHC2_251")), aes(ObsH_BP, ObsFst), shape = 0, size = 3) +
  geom_point(data = filter(OLinputRK, FstPval < 0.01), aes(ObsH_BP, ObsFst), shape = 2, size = 3) +
  theme_bw() +
  ggtitle("Red and Karluk") +
  scale_y_continuous(name = "Fst", breaks = seq(0,0.4,0.1)) + 
  scale_x_continuous(name = "Ho/(1-Fst)", breaks = seq(0,0.7,0.1)) +
  geom_line(aes(Het_BP, `0.01`), linetype = "dashed", data = OLSimInputRK) + 
  geom_line(aes(Het_BP, `0.99`), linetype = "dashed", data = OLSimInputRK) + 
  coord_cartesian(xlim = c(0, 0.7), ylim = c(-0.05,0.40)) +
  geom_text(data = filter(OLinputRK, FstPval < 0.01 | LocusName == c("MHC2_190", "MHC2_251")), 
            aes(ObsH_BP, ObsFst, label = LocusName), nudge_x = -0.035, nudge_y = 0.015, size = 4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin=unit(c(0.5,10,0.5,1.2),"cm"), 
        axis.title = element_text(size=14), axis.text = element_text(size=12), 
        plot.title = element_text(size = 14, vjust = -10, hjust = 0.1))


#plot Hs and Ar--------------------------------------------------------
Msat_Hs_plot <- ggplot(inputDivTest) + 
  geom_bar(aes(Lake, Msat_Hs, fill = Lake), color = "black", position = 'dodge',  stat = 'identity') + 
  scale_fill_manual(values = c("black", "gray", "white")) +
  theme_bw() + 
  ggtitle("Msat_Hs") +
  scale_x_discrete(limits = c("Karluk", "Red", "Frazer")) +
  scale_y_continuous(name = "Hs", breaks = seq(0,1,0.2)) +
  coord_cartesian(ylim = c(0,1)) +
  geom_text(aes(Lake, Msat_Hs, label=sprintf('%.3f',Msat_Hs)), vjust=-0.25, size = 2.5) +
  #annotate("text", x = 1, y = 2.25, label = c("P(Karluk>Frazer) = 0.001"), size = 4) + not needed becasue no significant test
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin=unit(c(0,0,0,0.5),"cm"), 
        axis.title = element_text(size=10), axis.text = element_text(size=8), 
        plot.title = element_text(size=12, vjust = 0, hjust = 0.05), legend.position = "none")

Msat_Ar_plot <- ggplot(inputDivTest) + 
  geom_bar(aes(Lake, Msat_Ar, fill = Lake), color = "black", position = 'dodge',  stat = 'identity') + 
  scale_fill_manual(values = c("black", "gray", "white")) +
  theme_bw() + 
  ggtitle("Msat_Ar") +
  scale_x_discrete(limits = c("Karluk", "Red", "Frazer")) +
  scale_y_continuous(name = "Ar", breaks = seq(0,18,2)) +
  coord_cartesian(ylim = c(0,18)) +
  geom_text(aes(Lake, Msat_Ar, label=Msat_Ar), vjust=-0.25, size = 2.5) +
  annotate("text", x = 1, y = 17.5, label = c("P(Karluk<=Frazer) = 0.001"), size = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin=unit(c(0,0,0,0.5),"cm"), 
        axis.title = element_text(size=10), axis.text = element_text(size=8), 
        plot.title = element_text(size=12, vjust = 0, hjust = 0.05), legend.position = "none")


#Comparing Hs and Ar-------------------------------------------------------------------------------
#used Fstat outside of R to compare Hs and Ar between Frazer and donor groups (Karluk beach, Red tributary)  
#convert genepop input file to Fstat input file - 
conversion(locinfile, "Fstat", "Kodiak_Nuc_Fstat.dat")

#Plot histogram of Hs and Ar from FSTAT output for SNPs and Msats in Karluk, Red, Frazer lakes.
inputDivTest <- read_csv("KodiakSockeye_DiversityTest3.csv") #read input as a tibble.
Hs_plot <- ggplot(inputDivTest) + 
  geom_bar(aes(Lake, Hs, fill = Lake), color = "black", position = 'dodge',  stat = 'identity') + 
  scale_fill_manual(values = c("black", "gray", "white")) +
  theme_bw() + 
  ggtitle("(a)") +
  scale_y_continuous(name = "Hs", breaks = seq(0,0.75,0.25)) +
  coord_cartesian(ylim = c(0,1)) +
  geom_text(aes(Lake, Hs, label=sprintf('%.3f',Hs)), vjust=-0.25, size = 2.5) +
  geom_text(aes(Lake, Hs, label = P_Hs), vjust=-1.5, size = 2, na.rm = T) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin=unit(c(0,0,0,0.5),"cm"), 
        axis.title = element_text(size=10), axis.text = element_text(size=8), 
        plot.title = element_text(size=12, vjust = 0, hjust = 0.05), legend.position = "none") +
  facet_grid(Marker ~ Habitat, scales = "free")

Ar_plot <- ggplot(inputDivTest) + 
  geom_bar(aes(Lake, Ar, fill = Lake), color = "black", position = 'dodge',  stat = 'identity') + 
  scale_fill_manual(values = c("black", "gray", "white")) +
  theme_bw() + 
  ggtitle("(b)") +
  scale_y_continuous(name = "Ar", breaks = seq(0,20,5)) +
  coord_cartesian(ylim = c(0,19)) +
  geom_text(aes(Lake, Ar, label=sprintf('%.2f',Ar)), vjust=-0.25, size = 2.5) +
  geom_text(aes(Lake, Ar, label = P_Ar), vjust=-1.5, size = 2, na.rm = T) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin=unit(c(0,0,0,0.75),"cm"), 
        axis.title = element_text(size=10), axis.text = element_text(size=8), 
        plot.title = element_text(size=12, vjust = 0, hjust = 0.05), legend.position = "none") +
  facet_grid(Marker ~ Habitat, scales = "free")

tiff('Figure3_FrazerPaper.tiff', height = 4.6, width = 3.4, units = 'in', compression = 'none', res = 600) #creates a high resolution TIFF file for publication
ggarrange(Hs_plot, Ar_plot, ncol = 1, nrow = 2)
dev.off()


#Test differentation at nuclear loci and mtDNA separately using Genepop package-----------------------------------------------
#Nuclear loci
test_diff(locinfile, genic = F, pairs = F, outputFile = "Kodiak_Gtest_genepop.txt", dememorization = 10000, batches = 100, iterations = 5000)
test_diff(locinfile, genic = F, pairs = T, outputFile = "Kodiak_GtestPW_genepop.txt", dememorization = 10000, batches = 100, iterations = 5000)
#mtDNA
test_diff(locinfile, genic = F, pairs = F, outputFile = "Kodiak_Gtest_mtDNA_genepop.txt", dememorization = 10000, batches = 100, iterations = 5000)
test_diff(locinfile, genic = F, pairs = T, outputFile = "Kodiak_GtestPW_mtDNA_genepop.txt", dememorization = 10000, batches = 100, iterations = 5000)


#Prepare a heat map for G-test results from mtDNA data----------------------------------------------------------------------------
#Prepare mtDNA PW G-test results for heat map 
PW_Gtest <- read_table("Kodiak_GtestPW_mtDNA_genepop_sub.txt", skip = 5, col_names = F) %>% 
  select(X1, X3, X6) %>% 
  rename(Pop1 = X1, Pop2 = X3, PW_Pval = X6) %>% 
  mutate(PW_Pval = recode(PW_Pval, "Highly sign." = "0.000000")) %>% 
  mutate(PW_Pval = str_remove(PW_Pval, "<"))
PW_Gtest_PopsCode <- unique(c(PW_Gtest$Pop1, PW_Gtest$Pop2))
PW_Gtest_PopNum <- as.character(c(1:17,19:33))
PW_Gtest_PopsCodeAbbr <- substr(PW_Gtest_PopsCode, 1,3)
PW_Gtest_PopsCodeAbbr[5] <- "FOM"
PW_Gtest_PopsCodeAbbrNum <- paste(PW_Gtest_PopNum, PW_Gtest_PopsCodeAbbr)
PW_Gtest_PopsCodeNum1 <- setNames(PW_Gtest_PopNum, PW_Gtest_PopsCode)
PW_Gtest_PopsCodeNum2 <- setNames(PW_Gtest_PopsCodeAbbrNum, PW_Gtest_PopsCode)
PW_Gtest$Pop1 <- recode(PW_Gtest$Pop1, !!!PW_Gtest_PopsCodeNum1)
PW_Gtest$Pop2 <- recode(PW_Gtest$Pop2, !!!PW_Gtest_PopsCodeNum2)
PW_Gtest$PW_Pval <- as.numeric(PW_Gtest$PW_Pval)
TempTib <- tibble(Pop1 = PW_Gtest_PopNum, Pop2 = PW_Gtest_PopsCodeAbbrNum, PW_Pval = rep(1, 32))
PW_Gtest <- bind_rows(TempTib, PW_Gtest)
PW_Gtest <- mutate(PW_Gtest, sig = ifelse(PW_Pval < 0.05, "1", "0"))
PopOrder <- unique(PW_Gtest$Pop1)
PopOrder2 <- unique(PW_Gtest$Pop2)
PW_Gtest$Pop1 <- factor(PW_Gtest$Pop1, levels = PopOrder)
PW_Gtest$Pop2 <- factor(PW_Gtest$Pop2, levels = PopOrder2)

#Plot heatmap of mtDNA pairwise g-test results
mtDNA_PW_GtestPlot <- ggplot(PW_Gtest, aes(Pop1, Pop2, fill = sig)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = c("white", "gray"), breaks = c(0,1), labels = c("NS", " <0.05"))  +
  theme_minimal() + 
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5, hjust = 0.5), 
        axis.text.y = element_text(angle = 0, vjust = 0.5, size = 5, hjust = 0.5), 
        axis.ticks = element_line(colour = "black", size = 0.25), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.25), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        legend.position = c(0.9, 0.9), 
        legend.text = element_text(size = 5, hjust = 0.5), 
        legend.title = element_blank(),  
        legend.key.size = unit(0.2, "cm"),
        plot.margin = margin(0,0,0,0.5, "cm"), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_discrete(position = 'top') + 
  scale_y_discrete(limits = rev(levels(factor(PW_Gtest$Pop2)))) +
  annotate("text", x = 1:32, y = 32:1, label = "-", size = 2) + 
  geom_rect(aes(xmin = 1 - 0.5, xmax = 10 + 0.5, ymin = 32 + 0.5, ymax = 22 + 0.5),
            fill = "transparent", color = "black", size = 0.25, linetype = 2) +
  geom_rect(aes(xmin = 11 - 0.5, xmax = 25 + 0.5, ymin = 22 + 0.5, ymax = 7 + 0.5),
            fill = "transparent", color = "black", size = 0.25, linetype = 2) +
  geom_rect(aes(xmin = 26 - 0.5, xmax = 31 + 0.5, ymin = 7 + 0.5, ymax = 1 + 0.5),
            fill = "transparent", color = "black", size = 0.25, linetype = 2) + 
  annotate("text", x = 5 + 0.5, y = 31, label = "Frazer", size = 2) + 
  annotate("text", x = 19 + 0.5, y = 21, label = "Karluk", size = 2) + 
  annotate("text", x = 29 + 0.5, y = 6, label = "Red", size = 2) +
  annotate("text", x = 16, y = 31, label = "mtDNA", size = 2.5)


Kodiak_mtDNAAMOVA_plot <- ggplot(filter(Kodiak_AMOVA, Class == "mtDNA"), aes(x = Strategy, y = Fct)) +
  geom_bar(stat = "identity", fill = "gray", color = "black") + 
  theme_bw() +
  scale_y_continuous(limits = c(0,0.035)) +
  scale_x_discrete(name = "Grouping Strategy") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")


#Alternative dapc plot with centroids symbols same as individual symbols and numbered-----------------------
dapc1_plot <- ggplot() + 
  geom_point(data = dapc1DF, aes(x = LD1, y = LD2, color = Group, shape = Group)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(),
        axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), legend.key = element_rect(size = 0.5),
        legend.key.size = unit(0.01, "lines"), legend.position = c(0.95, 0.70), legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5)) +
  scale_color_manual(values = myCol, labels = GrpLab) +
  scale_shape_manual(values = c(rep(1,3), rep(2,4), rep(5,2), 2), labels = GrpLab) +
  scale_y_continuous(limits = c(-5,7)) +
  scale_x_continuous(limits = c(-6,6)) +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  guides(color = guide_legend(override.aes = list(size = 1.75))) +
  #stat_ellipse(type = "t", show.legend = F, size = 0.75) +
  #annotate('label', x = dapc1$grp.coord[,1], y = dapc1$grp.coord[,2], label = c(1:10), cex = 4, fill = myCol, color = "white") +
  geom_point(data = dapc1GroupTib, aes(x = LD1, y = LD2, color = LakeGroup, shape = LakeGroup), size = 3, stroke = 2, show.legend = FALSE) +
  geom_text_repel(data = dapc1GroupTib, aes(x = LD1, y = LD2, label = c(1:10)), size = 4)


#Convert to genepop input file---------------------------------------------------------------------------------------
#rename input file from above (e.g., nuclear or mtDNA)
input2GP <- input_NucMtDNA_GoL
#Optional: remove "One_' from prefix of each locus name (so can distinguish loci in output files from GenePop analysis)
names(input2GP) <- gsub("One_","",names(input2GP),fixed=T) 
#Creates blank dataframe with same number of columns as inputGP and no. columns - 1 rows.
NewGP <-  as.tibble(matrix(nrow = length(input2GP) - 1, ncol = length(input2GP))) 
#Assign the first column of NewGP the column names of inputGP column 2 through end.  
NewGP[,1] <-  names(input2GP[,2:length(input2GP)]) 
#Assign the column names of inputGP to NewGP
names(NewGP) <-  names(input2GP) 
#Assign the first row, column 1 of NewGP the string "Genepop Input'.
NewGP[1,1] <- 'Genepop Input' 
#Create LociCount vector with one element that is the number of columns of NewGP minus 2.
LociCount <-  length(NewGP)-2 
input2GP[,-c(1,2)] <- mutate_if(input2GP[,-c(1,2)], is.character, str_replace_all, pattern = "/", replacement = "")
#Add a comma to each character in the column IndList in the dataframe inputGP.  The comma in front of the individual ID is needed for genepop. 
input2GP$Ind <- paste(input2GP$Ind, ',', sep = '') 
#Create a one row empty dataframe 
Pop <-  as.data.frame(t(c(rep(NA, length(input2GP))))) 
#Assign the column names of input2GP to Pop
names(Pop) <-  names(input2GP) 
#Assign the first row, column 1 of Pop the string 'Pop'
Pop[1] <-  'Pop' 
#Create a PopN vector with the names of each population from the the column input2GP$Pop in input2GP.
PopN <-  unique(input2GP$Pop)
#for loop appends the datafram NewGP with the Pop dataframe and each population, in turn, from the dataframe inputGP.
for(i in 1:length(PopN)) {
  tempGP <- filter(input2GP, Pop == PopN[i])
  NewGP <-  bind_rows(NewGP, Pop, tempGP)  
}
#Remove the column titled PopList from NewGP
NewGP <-  select(NewGP, -c(Pop)) 
#Write NewGP to file.
write.table(NewGP,"Kodiak_NucmtDNA_GP.txt",quote=FALSE,row.names=F,col.names=F,sep=" ",na="") 
write.table(PopN, "Kodiak_NucmtDNA_GP_PopN.txt", quote=FALSE,row.names=F,col.names=F,sep=" ",na="")
#remove temp variables used for loop
rm("i", "input2GP", "NewGP", "LociCount", "Pop", "PopN", "tempGP")


#Plot of One_U1004_183 allele frequency for each seasonal spawning time ecotype by lake------------------------------------
U1004_plot <- ggplot(U1004count, aes(x = PopGroup, y = freq)) + 
  geom_point(stat = "identity", size = 3) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 8)) +
  scale_color_manual(values = c("black")) +
  scale_shape_manual(values = c(19)) +
  scale_y_continuous(name = "Freq", limits = c(0,1)) +
  scale_x_discrete(labels = c("Frazer\nBeach", "Frazer\nOutlet", "Frazer\nTrib", 
                              "Karluk\nBeach", "Karluk\nOutlet", "Karluk\nTrib", "Karluk\nTrib", "Red\nBeach", "Red\nTrib")) + 
  annotate("text", x = U1004count$PopGroup, y = U1004count$freq + 0.07, label = c(rep("Mid", 3), rep("Late", 2), "Early", rep("Late", 2), "Early"), size = 2) +
  annotate("text", x = 2.5, y = 1.0, label = "One_U1004_183 Allele 001", size = 2.5)


# Pairwise Cavalli-Sforza and Edwards dist and G'st among all collections---------------------------------------------------------------------

#Use HierFstat to compute pairwise Cavalli-Sforza and Edwards genetic distance
Dchpw <- genet.dist(HF_NucMtDNA_GoL,diploid=TRUE,method="Dch")
#convert Dchpw to dist list and add columns for filtering
DchpwList <- dist2list(Dchpw) %>% 
  select(row, col, value) %>% 
  mutate("R" = as.numeric(substr(row,4,5)), "C" = as.numeric(substr(col,4,5)))
#compute mean pairwise Dch for each Red Lake collection overall Frazer Lake collections
DchpwListSummarize <- DchpwList %>% 
  filter(R != C & !(R > C)) %>% 
  filter(R %in% c(1:10), C %in% c(27:32)) %>% 
  group_by(C) %>% summarize(mean(value))

#Use mmod to estimate pairwise G'st
HedGSTpw <- pairwise_Gst_Hedrick(AD_NucMtDNA_GoL, linearized = FALSE)
#Use mmod to estimate G'st overall and for each locus
HedGST <- Gst_Hedrick(AD_NucMtDNA_GoL)
#Gst by locus
GstByLocus <- as_tibble(HedGST$per.locus) %>% 
  mutate(Locus = LociList_NucMtDNA_GoL) %>% 
  rename(Gst = value) %>% 
  select(Locus, Gst)


# Multi-dimensional scaling of PW Cavalli  estimates for all collections-------------------------------
MDSDch <- cmdscale(Dchpw,eig=FALSE,k=3) %>% 
  as_tibble() %>% 
  rename(Dimension_1 = V1, Dimension_2 = V2, Dimension_3 = V3) %>% 
  mutate(CollectionNo = CollectionNo, Ecotype = Ecotype) %>% 
  select(CollectionNo, Ecotype, everything())

#convert WCpw to dist matrix--------------
WCpwDist <- as.dist(WCpw)

#convert WCpwDist to dist list and add columns for filtering------------
WCpwList <- dist2list(WCpwDist) %>% 
  select(row, col, value) %>% 
  mutate("R" = as.numeric(substr(row,4,5)), "C" = as.numeric(substr(col,4,5))) %>% 
  filter(!(R > C))

#compute mean pairwise Fst for each Red Lake collection overall Frazer Lake collections---------------------
WCpwListSummarize <- WCpwList %>% 
  filter(R != C & !(R > C)) %>% 
  filter(R %in% c(1:10), C %in% c(27:32)) %>% 
  group_by(C) %>% summarize(mean(value))


#Scatter plot of MDS of PW Cavalli-Sforza and Edwards distance estimates--------------------------
PWDchMDSplot <- ggplot(MDSDch, aes(x = Dimension_1, y = Dimension_2, color = Ecotype, shape = Ecotype)) + 
  #Data visualization  
  geom_point(stroke = 1) +
  #layout
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=10), axis.text = element_text(size=8), 
        legend.key = element_rect(size = 0.5), legend.key.size = unit(0.1, "lines"), 
        legend.position = c(0.95, 0.80), legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5), legend.spacing.y = unit(0.25, 'lines'), 
        legend.background = element_rect(fill='transparent')) +
  #overrides
  scale_color_manual(values = myColEcotype) +
  scale_shape_manual(values = c(rep(19,3), rep(2,4), rep(5,3), 25)) +
  scale_y_continuous(limits = c(-0.02,0.03), breaks = seq(-0.02,0.03,0.01), name = "Dimension 2") +
  scale_x_continuous(limits = c(-0.02,0.03), breaks = seq(-0.02,0.03,0.01), name = "Dimension 1") +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  guides(color = guide_legend(override.aes = list(size = 1.5), byrow = TRUE)) +
  annotate('text', x = 0.015, y = -0.007, label = "Karluk, Late\nKBL,KOL,KTL", size = 2,  color = "black") +
  annotate('text', x = 0.015, y = 0.005, label = "Karluk, Early\nKTE", size = 2,  color = "black") +
  annotate('text', x = -0.010, y = -0.006, label = "Red, Late\nRBL", size = 2,  color = "black") +
  annotate('text', x = -0.019, y = 0.007, label = "Red, Early\nRTE", size = 2,  color = "black") +
  annotate('text', x = -0.007, y = 0.005, label = "Frazer\nFBM,FOM,FTM", size = 2,  color = "black") +
  annotate('text', x = 0.021, y = 0.028, label = "RUL", size = 2,  color = "black")  +
  annotate('text', x = MDSDch$Dimension_1[13], y = MDSDch$Dimension_2[13] - 0.002, label = "13", size = 2,  color = "black") +
  annotate('text', x = MDSDch$Dimension_1[31], y = MDSDch$Dimension_2[31] + 0.002, label = "31", size = 2,  color = "black") +
  annotate('text', x = MDSDch$Dimension_1[32], y = MDSDch$Dimension_2[32] + 0.002, label = "33", size = 2,  color = "black") +
  annotate('text', x = -0.018, y = 0.03, label = expression("D"[CH]), size = 4,  color = "black")


MDS_PWFstCSE <- ggarrange(PWFstMDSplot, PWCSEMDSplot, ncol = 1, nrow = 2, 
                          heights = c(8, 8, 1), common.legend = TRUE, legend = "bottom")


#Optional Use poppr to compute pairwise Edwards dist----------------------------------------
EdwardsPW <- edwards.dist(Gen_NucMtDNA_GoL)

#Optional Multi-dimensional scaling of PW Edwards distance estimates for all collections-------------------------------
#Perform MDS and plot results
MDSEdwards<-cmdscale(EdwardsPW,eig=FALSE,k=3) %>% 
  as_tibble(.name_repair = NULL) %>% 
  rename(Dimension_1 = V1, Dimension_2 = V2, Dimension_3 = V3) %>% 
  mutate(SampNo = SampNo, Ecotype = Ecotype) %>% 
  select(SampNo, Ecotype, everything())

#Optional Scatter plot of MDS of PW Edwards distance estimates--------------------------
PWEdwardsMDSplot <- ggplot(MDSEdwards, aes(x = Dimension_1, y = Dimension_2, color = Ecotype, shape = Ecotype)) + 
  #Data visualization  
  geom_point(stroke = 1) +
  #layout
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=10), axis.text = element_text(size=8), 
        legend.key = element_rect(size = 0.5), legend.key.size = unit(0.1, "lines"), 
        legend.position = c(0.95, 0.80), legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5), legend.spacing.y = unit(0.25, 'lines'), 
        legend.background = element_rect(fill='transparent')) +
  #overrides
  scale_color_manual(values = myColEcotype) +
  scale_shape_manual(values = c(rep(19,3), rep(2,4), rep(5,3), 25)) +
  scale_y_continuous(limits = c(-0.05,0.11), breaks = seq(-0.05,0.11,0.02), name = "Dimension 2") +
  scale_x_continuous(limits = c(-0.07,0.07), breaks = seq(-0.07,0.07,0.02), name = "Dimension 1") +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  guides(color = guide_legend(override.aes = list(size = 1.5), byrow = TRUE)) +
  annotate('text', x = 0.063, y = -0.027, label = "Karluk, Late\nKBL,KOL,KTL", size = 2,  color = "black") +
  annotate('text', x = 0.063, y = 0.033, label = "Karluk, Early\nKTE", size = 2,  color = "black") +
  annotate('text', x = -0.050, y = -0.032, label = "Red, Late\nRBL", size = 2,  color = "black") +
  annotate('text', x = -0.050, y = 0.035, label = "Red, Early\nRTE", size = 2,  color = "black") +
  annotate('text', x = -0.063, y = -0.008, label = "Frazer\nFBM,FOM,FTM", size = 2,  color = "black") +
  annotate('text', x = 0.040, y = 0.105, label = "RUL", size = 2,  color = "black")  +
  annotate('text', x = MDSEdwards$Dimension_1[13], y = MDSEdwards$Dimension_2[13] - 0.004, label = "13", size = 2,  color = "black") +
  annotate('text', x = MDSEdwards$Dimension_1[31], y = MDSEdwards$Dimension_2[31] + 0.006, label = "31", size = 2,  color = "black") +
  annotate('text', x = MDSEdwards$Dimension_1[32], y = MDSEdwards$Dimension_2[32] + 0.006, label = "33", size = 2,  color = "black") +
  annotate('text', x = -0.065, y = 0.11, label = "CSE", size = 4,  color = "black")

#Optional Scatter plot of MDS of PW Fst estimates for Frazer and donor populations only--------------------------
#Color for each ecotype in plot
myColEcotype = c("dodgerblue4", 'dodgerblue', 'deepskyblue', 'springgreen3',"gold2", "slateblue")

#Scatter plot of MDS of PW Fst estimates including locus One_U1004_183
Figure4x <- ggplot(MDSFst, aes(x = Dimension_1, y = Dimension_2, color = Ecotype, shape = Ecotype)) + 
  #Data visualization  
  geom_point(stroke = 1) +
  #layout
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=8), 
        legend.key = element_rect(linewidth = 0.5), legend.key.size = unit(0.1, "lines"), 
        legend.position = "bottom", legend.title = element_text(size = 7), 
        legend.text = element_text(size = 6), legend.spacing.y = unit(0.25, 'lines'), 
        legend.background = element_rect(fill='transparent')) +
  #overrides
  scale_color_manual(values = myColEcotype) +
  #  scale_shape_manual(values = c(rep(19,3), rep(2,4), rep(5,3), 25)) +
  scale_shape_manual(values = c(17, 15, 18, 2, 2, 5)) +
  scale_y_continuous(limits = c(-0.04,0.05), breaks = seq(-0.04,0.05,0.01), name = "Dimension 2") +
  scale_x_continuous(limits = c(-0.03,0.07), breaks = seq(-0.03,0.07,0.01), name = "Dimension 1") +
  #  scale_y_continuous(limits = c(-0.02,0.02), breaks = seq(-0.02,0.02,0.01), name = "Dimension 2") +
  #  scale_x_continuous(limits = c(-0.02,0.03), breaks = seq(-0.02,0.03,0.01), name = "Dimension 1") +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  guides(color = guide_legend(override.aes = list(size = 1.5), nrow = 1)) +
  #annotate('text', x = 0.019, y = 0.022, label = "Karluk, Late\nKBL,KOL,KTL", size = 2,  color = "black") +
  #annotate('text', x = 0.022, y = 0.005, label = "Karluk, Early\nKTE", size = 2,  color = "black") +
  #annotate('text', x = -0.017, y = 0.005, label = "Red, Late\nRBL", size = 2,  color = "black") +
  #annotate('text', x = -0.021, y = -0.028, label = "Red, Early\nRTE", size = 2,  color = "black") +
  #annotate('text', x = -0.013, y = -0.014, label = "Frazer\nFBM,FOM,FTM", size = 2,  color = "black") +
  #annotate('text', x = 0.019, y = 0.022, label = "Karluk, Late", size = 2,  color = "black") +
  #annotate('text', x = 0.022, y = 0.005, label = "Karluk, Early", size = 2,  color = "black") +
  #annotate('text', x = -0.017, y = 0.005, label = "Red, Late", size = 2,  color = "black") +
  #annotate('text', x = -0.021, y = -0.028, label = "Red, Early", size = 2,  color = "black") +
  #annotate('text', x = -0.010, y = -0.011, label = "Frazer", size = 2,  color = "black") +
  #annotate('text', x = 0.051, y = -0.039, label = "Ruth", size = 2,  color = "black") +
annotate('text', x = MDSFst$Dimension_1[11], y = MDSFst$Dimension_2[11] + 0.003, label = "13", size = 2,  color = "black") +
  annotate('text', x = MDSFst$Dimension_1[12], y = MDSFst$Dimension_2[12] + 0.003, label = "31", size = 2,  color = "black") +
  annotate('text', x = MDSFst$Dimension_1[13], y = MDSFst$Dimension_2[13] + 0.003, label = "33", size = 2,  color = "black") +
  annotate('text', x = -0.027, y = 0.05, label = "A", size = 4,  color = "black", parse = TRUE)


#Phylogram of all collections using Adegenet, Poppr and Ape--------
#Use Adegenet file with all loci 
#Prepare for use in adegenet, hierfsat packages
#Vector of populations from nuclear data (for adegenet)
PopList_NucMtDNA_GoL <- input_NucMtDNA_GoL$Pop #List of populations for Adegent
#Vector of individuals from nuclear data (for adegenet)
IndList_NucMtDNA_GoL <- input_NucMtDNA_GoL$Ind #List of individuals from Adegent
#Vector of nuclear loci
LociList_NucMtDNA_GoL <- names(select(input_NucMtDNA_GoL, -c(Pop, Ind)))
#Use 'select' in dplyr package to remove columns Ind, Pop.
input_NucMtDNA_GoL4AD <- select(input_NucMtDNA_GoL, -c(Ind, Pop)) 
#Convert to adegenet input file
AD_NucMtDNA_GoL <- df2genind(input_NucMtDNA_GoL4AD,sep = '/', ncode=6, ind.names=IndList_NucMtDNA_GoL, pop=PopList_NucMtDNA_GoL, NA.char="000/000", type="codom") 
#Convert to hierfstat input file
HF_NucMtDNA_GoL <- genind2hierfstat(AD_NucMtDNA_GoL) 
#Use Adegenet to convert genind object to genpop object
Gen_NucMtDNA_GoL <- genind2genpop(AD_NucMtDNA_GoL)
#Use Poppr to make NJ tree based on Nei's distance and showing bootstrap values for nodes
NJtreeNei <- aboot(Gen_NucMtDNA_GoL, strata = NULL, tree = "nj", distance = "nei.dist", 
                   sample = 1000, cutoff = 75, showtree = TRUE, missing = "mean", mcutoff = 0, 
                   quiet = FALSE, root = FALSE)


#Alternative phylogram of all collections from each lake--------------------
# Define colors for branch tip labels - blue for Frazer Lake collections and black for others
NJPopCols <- c(rep("blue", 10), rep("black", 2), "green", rep("black", 16), "orange", "black", "gold4")
#Use Ape to plot unrooted tree from Poppr - used the plotbreaklongedges function instead of plot to break the long branch length for RUL33 (the argument n=1 refers to one branch length)
tiff('Figure3_FrazerPaper.tiff', height = 3.5, width = 3.5, units = 'in', compression = 'none', res = 600) #creates a high resolution TIFF file for publication
plotBreakLongEdges(NJtreeNei, n = 1, type = "unrooted", show.node.label = TRUE, edge.width = 0.5, 
                   tip.color = NJPopCols, cex = 0.4, font = 1, lab4ut = "axial", no.margin = T)
#this code can be used if not breaking long branch length
#plot(NJtreeNei, type = "unrooted", show.node.label = TRUE, edge.width = 0.5, 
#                   tip.color = NJPopCols, cex = 0.4, font = 1, lab4ut = "axial", no.margin = T)
add.scale.bar(cex = 0.4)
text(0,0.018,"Red and Frazer\nLakes", cex = 0.6, pos = 4)
text(0.02, 0.007, "Karluk Lake", cex = 0.6)
dev.off()

#test code
plotBreakLongEdges(testtree, n = 1, type = "unrooted", show.node.label = TRUE, edge.width = 0.5, cex = 0.3, font = 1, lab4ut = "axial", no.margin = T)
GenDistMatrix <- dist.genpop(Gen_NucMtDNA_GoL,1)
NJ_GDM <- nj(GenDistMatrix)

# unused figure for within-lake AMOVA
Kodiak_AMOVA_InLakes_plot <- ggplot(Kodiak_AMOVA_InLakes) + 
  geom_bar(aes(Fstatistic, Val), fill = "gray", color = "black", position = 'dodge',  stat = 'identity') + 
  #scale_fill_manual(values = c("black", "gray", "white")) +
  theme_bw() + 
  #ggtitle("(a)") +
  scale_y_continuous(name = "Val", breaks = seq(0,0.02,0.01)) +
  coord_cartesian(ylim = c(0,0.030)) +
  scale_x_discrete(limits = c("Fst", "Fsc", "Fct")) +
  geom_text(aes(Fstatistic, Val, label=sprintf('%.3f',Val)), vjust=-0.25, size = 2.5) +
  geom_text(aes(Fstatistic, Val, label = ifelse(Signif == "ns", "", Signif)), vjust = -1.0, hjust = -3.5, size = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=10), axis.text = element_text(size=8), 
        plot.title = element_text(size=12, vjust = 0, hjust = 0.05), legend.position = "none") +
  facet_grid(factor(Lake, levels = c("Karluk", "Red", "Frazer")) ~ U1004_YN, scales = "free")
