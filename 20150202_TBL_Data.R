## Compare Biggest Loser Genotype Frequencies w/ 1KG Data ##
## August 25, 2014 ##
## Kristopher Standish ##

############################################################################################################
## OPEN UP ALL THE FILES ###################################################################################
############################################################################################################

start <- proc.time()
library(xlsx)
library(glmnet)

######################################################
## SET UP PATHS ######################################
######################################################

## Path To Save Plots/Results to
DATE <- "20150202"
PathToOut <- paste("/Users/kstandis/Data/TBL/Plots/",DATE,"_Plots",sep="")

## Path to Old Data (8/25)
PathToOldData <- "/Users/kstandis/Data/TBL/Data/FILENAME"

## Path to New Set of Data
PathToNewData <- "/Users/kstandis/Data/TBL/Data/20141001/FILENAME"

## Path to New Set of 1KG VCF files
PathTo1KG <- "/Users/kstandis/Data/TBL/Data/20141001/1KG/FILENAME"

######################################################
## LOAD DATA #########################################
######################################################

## Load Old Genotype Data (8/25)
AFR_825 <- read.table(gsub("FILENAME","1KG_AFR_Tab.txt",PathToOldData), sep="\t",header=T, colClasses=c("numeric","numeric","character","factor","factor",rep("numeric",3)))
EUR_825 <- read.table(gsub("FILENAME","1KG_EUR_Tab.txt",PathToOldData), sep="\t",header=T, colClasses=c("numeric","numeric","character","factor","factor",rep("numeric",3)))
ASN_825 <- read.table(gsub("FILENAME","1KG_ASN_Tab.txt",PathToOldData), sep="\t",header=T, colClasses=c("numeric","numeric","character","factor","factor",rep("numeric",3)))
SAN_825 <- read.table(gsub("FILENAME","1KG_SAN_Tab.txt",PathToOldData), sep="\t",header=T, colClasses=c("numeric","numeric","character","factor","factor",rep("numeric",3)))
AMR_825 <- read.table(gsub("FILENAME","1KG_AMR_Tab.txt",PathToOldData), sep="\t",header=T, colClasses=c("numeric","numeric","character","factor","factor",rep("numeric",3)))
ALL_825 <- read.table(gsub("FILENAME","1KG_ALL_Tab.txt",PathToOldData), sep="\t",header=T, colClasses=c("numeric","numeric","character","factor","factor",rep("numeric",3)))
TBL_825 <- read.table(gsub("FILENAME","1KG_TBL_Tab_Reformat2.txt",PathToOldData), sep="\t",header=T)
start-proc.time()

## Load rsID/Gene/Phenotype Info (10/01)
PGR_DIET <- read.xlsx(gsub("FILENAME","FIT genes.xlsx",PathToNewData), sheetName="Diet", colIndex=3:6, rowIndex=2:69, as.data.frame=T, header=T)
PGR_NUTR <- read.xlsx(gsub("FILENAME","FIT genes.xlsx",PathToNewData), sheetName="Nutrition", colIndex=3:6, rowIndex=2:11, as.data.frame=T, header=T)
PGR_EXER <- read.xlsx(gsub("FILENAME","FIT genes.xlsx",PathToNewData), sheetName="Exercise", colIndex=3:6, rowIndex=3:14, as.data.frame=T, header=T)
PGR_META <- read.xlsx(gsub("FILENAME","FIT genes.xlsx",PathToNewData), sheetName="Metabolic Health", colIndex=3:6, rowIndex=2:55, as.data.frame=T, header=T)
start-proc.time()

## Load Phenotype/Genotype/Clinical Data for Patients (10/01)
PHENO <- read.table(gsub("FILENAME","Phenotypic-Outcomes-Format-2.txt",PathToNewData), header=T, sep="\t" ) ; start-proc.time()
GENO <- read.table(gsub("FILENAME","Genotypes-Format-2.txt",PathToNewData), header=T, sep="\t", fill=T, comment.char="" ) ; start-proc.time()
CLIN.l <- read.xlsx(gsub("FILENAME","BL all sample list_130426_FINAL_pg1.xlsx",PathToNewData), sheetIndex=1, colIndex=1:139, rowIndex=1:289, as.data.frame=T, header=T)
CLIN_KEY <- data.frame( NM_1=colnames(CLIN.l), NM_2=c(CLIN.l[1,],recursive=T), row.names=NULL )
CLIN <- CLIN.l[2:nrow(CLIN.l),]
POOL <- read.table(gsub("FILENAME","SNP_Pheno_Pool_Summary.txt",PathToNewData), header=T, sep="\t", fill=T )
start-proc.time()

## Load 1KG files (10/01) 
 # Pedigree Info
PED <- read.table(gsub("FILENAME","1KG_Pedigree_20130606_g1k.ped",PathTo1KG), header=T, sep="\t")
INFO_PED <- read.xlsx(gsub("FILENAME","1KG_20130606_sample_info.xlsx",PathTo1KG), sheetName="Sample Info", colIndex=1:15, rowIndex=1:3501, as.data.frame=T, header=T)
start-proc.time()
 # VCF Files
VCF_ALL <- read.table(gsub("FILENAME","1KG_All_Vars.vcf",PathTo1KG), sep="\t", header=T, skip=254,comment.char="")
VCF_AFR <- read.table(gsub("FILENAME","1KG_AFR.recode.vcf",PathTo1KG), sep="\t", header=T, skip=254,comment.char="")
VCF_AMR <- read.table(gsub("FILENAME","1KG_AMR.recode.vcf",PathTo1KG), sep="\t", header=T, skip=254,comment.char="")
VCF_ASN <- read.table(gsub("FILENAME","1KG_ASN.recode.vcf",PathTo1KG), sep="\t", header=T, skip=254,comment.char="")
VCF_SAN <- read.table(gsub("FILENAME","1KG_SAN.recode.vcf",PathTo1KG), sep="\t", header=T, skip=254,comment.char="")
VCF_EUR <- read.table(gsub("FILENAME","1KG_EUR.recode.vcf",PathTo1KG), sep="\t", header=T, skip=254,comment.char="")
start-proc.time()
 # Ancestry Key
KEY <- read.table(gsub("FILENAME","Panel_Key.txt",PathTo1KG), header=T)

## Load

############################################################################################################
## FILTER TBL DATA #########################################################################################
############################################################################################################

######################################################
## FILTER THRU FAMILY RELATIONS (CLIN) ###############
######################################################

## Which relations do we NOT know?
RM_UNKNOWN <- grep( "unknown", CLIN$family, ignore.case=T )
RM_BLFILL <- grep( "BL fill in", CLIN$family, ignore.case=T )
RM_QUEST <- grep( "?", CLIN$family, fixed=T )

RM_DUNNO <- Reduce( union, list(RM_UNKNOWN, RM_BLFILL, RM_QUEST) )

## Remove Second Generation
 # Remove all Sons & Daughters
RM_SONS <- grep( "son", CLIN$family, ignore.case=T )
RM_DAU <- grep( "daug", CLIN$family, ignore.case=T )
 # Remove Cousins
RM_COUS <- grep( "cous", CLIN$family, ignore.case=T )

RM_KIDS <- Reduce( union, list(RM_SONS, RM_DAU, RM_COUS) )

## Mark Brothers, Sisters, and Twins
MK_TWINS <- grep( "twin", CLIN$family, ignore.case=T ) # 1 set of Twins
MK_BRO <- grep( "broth", CLIN$family, ignore.case=T ) # 1 pair of Brothers
MK_SIS <- grep( "sist", CLIN$family, ignore.case=T ) # 0 sisters

MK_SIBS <- Reduce( union, list(MK_TWINS, MK_BRO, MK_SIS) )

RM_SIBS <- c( sample(MK_TWINS,1), sample(MK_BRO,1) ) # Based on 1 set of twins, 1 pair of brothers, 0 sisters
KP_SIBS <- setdiff( MK_SIBS, RM_SIBS )

## Keep Top Generation
 # Keep Fathers & Mothers
KP_DAD <- union( grep( "dad", CLIN$family, ignore.case=T ), grep( "fath", CLIN$family, ignore.case=T ) )
KP_MOM <- union( grep( "mom", CLIN$family, ignore.case=T ), grep( "moth", CLIN$family, ignore.case=T ) )

 # Keep Husbands & Wives (& Spouses)
KP_HUSB <- grep( "husb", CLIN$family, ignore.case=T )
KP_WIFE <- grep( "wife", CLIN$family, ignore.case=T )
KP_SPOUSE <- grep( "spouse", CLIN$family, ignore.case=T )

KP_PAR <- Reduce( union, list(KP_DAD, KP_MOM, KP_HUSB, KP_WIFE, KP_SPOUSE) )

## Keep people Known to have NO relations?
KP_NO <- grep( "no", CLIN$family, ignore.case=T ) # This actually includes a couple "unkNOwn"...will remove later

## Compile which samples to Keep and Remove
KP.1 <- Reduce( union, list(KP_SIBS, KP_PAR, KP_NO) )
RM <- Reduce( union, list(RM_DUNNO, RM_KIDS, RM_SIBS) )

KP <- setdiff( KP.1, RM )

## Any samples NOT accounted for?
MISSING <- setdiff(1:nrow(CLIN),union(KP,RM)) # Should be good...

######################################################
## FILTER THRU SELF-REPORTED ETHNICITY (CLIN) ########
######################################################

## Check out all the different Ethnicities in the cohort
table( CLIN$Ethnicity )

## Get all Unknown Ethnicities
ETH_DUNNO <- grep( "BL fill in", CLIN$Ethnicity, ignore.case=T )

## Get all "Caucasians"
ETH_CAUC.ALL <- grep( "aucasian", CLIN$Ethnicity, ignore.case=T )
 # Get all Caucasians + Other
ETH_CAUC.PLUS1 <- grep( "aucasian-", CLIN$Ethnicity, fixed=T )
ETH_CAUC.PLUS2 <- grep( "/caucasian", CLIN$Ethnicity, fixed=T )
 # Pull out ONLY full Caucasians
ETH_CAUC.ONLY <- setdiff( ETH_CAUC.ALL, union(ETH_CAUC.PLUS1,ETH_CAUC.PLUS2) )

## Get Hispanic
ETH_HISP <- grep( "hispanic", CLIN$Ethnicity, ignore.case=T )

## Get African/African-American
ETH_AFR <- grep( "african", CLIN$Ethnicity, ignore.case=T )

## Get "Other"
ETH_OTH <- setdiff( 1:nrow(CLIN), Reduce(union,list(ETH_CAUC.ONLY,ETH_HISP,ETH_AFR,ETH_DUNNO)) )

ETH <- list( Dunno=ETH_DUNNO, Cauc=ETH_CAUC.ONLY, Hisp=ETH_HISP, Afr=ETH_AFR, Other=ETH_OTH )

######################################################
## GET RID OF MOSTLY MISSING GENOTYPES ###############
######################################################

## Identify Columns in GENO file that are largely missing
SNPS_TO_REMOVE <- c()
for ( c in 2:ncol(GENO) ) {
	if ( length(which( GENO[,c]=="-/-")) > 0.9*nrow(GENO) ) {
		SNPS_TO_REMOVE <- c( SNPS_TO_REMOVE, c )
	}
}

## Remove Mostly Missing SNPs from Table
GENO.R <- GENO[,-SNPS_TO_REMOVE]

######################################################
## GET FINAL TBL SAMPLE LISTS ########################
######################################################

## Keep + Caucasians
KP_CAUC <- intersect( KP, ETH$Cauc )
KP_CAUC_SAMPS <- CLIN$Accession.number[KP_CAUC]

## Get Intersection of GENO, PHENO, CLIN sample names
KP_INT <- Reduce( intersect, list(KP_CAUC_SAMPS,PHENO$Accession_No,GENO$Accession..) )

## Pull out Sample Rows for PHENO/GENO/CLIN Tables
PHENO.C <- PHENO[ which(PHENO$Accession_No %in% KP_INT) , ]
GENO.C <- GENO.R[ which(GENO.R$Accession.. %in% KP_INT) , ]
# GENO.C <- GENO[ which(GENO$Accession.. %in% KP_INT) , ]
CLIN.C <- CLIN[ which(CLIN$Accession.number %in% KP_INT) , ]

############################################################################################################
## FILTER 1KG DATA #########################################################################################
############################################################################################################

######################################################
## GET UNRELATED SAMPLES #############################
######################################################

## Take only "mother", "father", and "unrel"/"unrels" (& NA)
 # PED File
PED_MOM <- which( PED$Relationship=="mother" )
PED_DAD <- which( PED$Relationship=="father" )
PED_UNR <- which( PED$Relationship=="unrels" | PED$Relationship=="unrel" )
PED_KP <- Reduce( union, list(PED_MOM,PED_DAD,PED_UNR) )
 # INFO_PED File
IPED_MOM <- which( INFO_PED$Relationship=="mother" )
IPED_DAD <- which( INFO_PED$Relationship=="father" )
IPED_UNR <- which( INFO_PED$Relationship=="unrels" | INFO_PED$Relationship=="unrel" )
IPED_NA <- which(is.na( INFO_PED$Relationship ))
IPED_KP <- Reduce( union, list(IPED_MOM,IPED_DAD,IPED_UNR,IPED_NA) )

## Specify Sample Names for PED_KP/IPED_KP
PED_KP_SAMPS <- PED$Individual.ID[PED_KP]
IPED_KP_SAMPS <- INFO_PED$Sample[IPED_KP]
KP_SAMPS <- intersect( IPED_KP_SAMPS,PED_KP_SAMPS )

######################################################
## SANITY CHECK ######################################
######################################################

## Sanity Check #1
table( table( PED$Family.ID[ which(PED$Individual.ID %in% KP_SAMPS) ] ) )
names( which( table( PED$Family.ID[ which(PED$Individual.ID %in% KP_SAMPS) ] ) > 2 ) )
## Sanity Check #2
length( which( PED$Siblings[ which(PED$Individual.ID %in% KP_SAMPS) ] != 0 ) )
head( PED[ intersect( which(PED$Individual.ID %in% KP_SAMPS), which( PED$Siblings!=0 ) ) , ] )

## Loop through each KP_SAMPS sample and check each familial relation
REMOVE_RELATED <- function(KP_SAMPS) {
	HAVE_FAM_PED <- HAVE_FAM_PED_RAND <- array(,c(0,6)) ; colnames(HAVE_FAM_PED) <- colnames(PED)[c(2,3:4,9:11)]
	HAVE_FAM_IPED <- HAVE_FAM_IPED_RAND <- array(,c(0,9)) ; colnames(HAVE_FAM_IPED) <- colnames(INFO_PED)[c(1,7:14)]
	PROBLEM_PED <- c() ; PROBLEM_PED_RAND <- c()
	PROBLEM_IPED <- c() ; PROBLEM_IPED_RAND <- c()
	for ( samp in sample(KP_SAMPS,length(KP_SAMPS),replace=F) ) {
		## PED Table ##
		PED_ROW <- which( PED$Individual.ID == samp )
		if ( any( PED[PED_ROW,c(3:4,9:11)]!=0 ) & !(samp %in% PROBLEM_PED) ) {
			## Compile Info about person w/ family
			HAVE_FAM_PED <- rbind(HAVE_FAM_PED, PED[PED_ROW,c(2,3:4,9:11)])
			## Check if any family members are in "KP_SAMPS"
			# Split strings
			FAM_COMPILE <- c()
			for ( column in c(3:4,9:11) ) {
				if ( PED[PED_ROW,column]!=0 ) {
					Strip <- gsub(" ","",as.character(PED[PED_ROW,column]))
					FAM_COMPILE <- c( FAM_COMPILE, strsplit( Strip, "," )[[1]] )
				}
			}
			WHICH_IN <- KP_SAMPS[which(KP_SAMPS %in% FAM_COMPILE)]
			PROBLEM_PED <- c( PROBLEM_PED, WHICH_IN )
		} # Close PED loop
		## PED Table (Rand) ##
		if ( any( PED[PED_ROW,c(3:4,9:11)]!=0 ) & !(samp %in% PROBLEM_PED_RAND) ) {
			## Compile Info about person w/ family
			HAVE_FAM_PED_RAND <- rbind(HAVE_FAM_PED_RAND, PED[PED_ROW,c(2,3:4,9:11)])
			## Check if any family members are in "KP_SAMPS"
			# Split strings
			FAM_COMPILE <- c()
			for ( column in c(3:4,9:11) ) {
				if ( PED[PED_ROW,column]!=0 ) {
					Strip <- gsub(" ","",as.character(PED[PED_ROW,column]))
					FAM_COMPILE <- c( FAM_COMPILE, strsplit( Strip, "," )[[1]] )
				}
			}
			WHICH_IN_RAND <- c( KP_SAMPS[which(KP_SAMPS %in% FAM_COMPILE)], samp )
			PROBLEM_PED_RAND <- c( PROBLEM_PED_RAND, sample(WHICH_IN_RAND,length(WHICH_IN_RAND)-1) )
		} # Close PED (Rand) loop

		## INFO_PED Table ##
		IPED_ROW <- which( INFO_PED$Sample == samp )
		if ( any( !is.na(INFO_PED[IPED_ROW,c(7:14)]) ) & !(samp %in% PROBLEM_IPED) ) { # HAVE_FAM_IPED <- c(HAVE_FAM_IPED, INFO_PED[IPED_ROW,c(1,7:14)]) }
			## Compile Info about person w/ family
			HAVE_FAM_IPED <- rbind(HAVE_FAM_IPED, INFO_PED[IPED_ROW,c(1,7:14)])
			## Check if any family members are in "KP_SAMPS"
			# Split strings
			FAM_COMPILE <- c()
			for ( column in 7:14 ) {
				if ( !is.na(INFO_PED[IPED_ROW,column]) ) {
					Strip <- gsub(" ","",as.character(INFO_PED[IPED_ROW,column]))
					FAM_COMPILE <- c( FAM_COMPILE, strsplit( Strip, "," )[[1]] )
				}
			}
			WHICH_IN <- KP_SAMPS[which(KP_SAMPS %in% FAM_COMPILE)]
			PROBLEM_IPED <- c( PROBLEM_IPED, WHICH_IN )
		} # Close INFO_PED loop
		## INFO_PED Table (Rand) ##
		if ( any( !is.na(INFO_PED[IPED_ROW,c(7:14)]) ) & !(samp %in% PROBLEM_IPED_RAND) ) { # HAVE_FAM_IPED <- c(HAVE_FAM_IPED, INFO_PED[IPED_ROW,c(1,7:14)]) }
			## Compile Info about person w/ family
			HAVE_FAM_IPED_RAND <- rbind(HAVE_FAM_IPED_RAND, INFO_PED[IPED_ROW,c(1,7:14)])
			## Check if any family members are in "KP_SAMPS"
			# Split strings
			FAM_COMPILE <- c()
			for ( column in 7:14 ) {
				if ( !is.na(INFO_PED[IPED_ROW,column]) ) {
					Strip <- gsub(" ","",as.character(INFO_PED[IPED_ROW,column]))
					FAM_COMPILE <- c( FAM_COMPILE, strsplit( Strip, "," )[[1]] )
				}
			}
			WHICH_IN_RAND <- c( KP_SAMPS[which(KP_SAMPS %in% FAM_COMPILE)], samp )
			PROBLEM_IPED_RAND <- c( PROBLEM_IPED_RAND, sample(WHICH_IN_RAND,length(WHICH_IN_RAND)-1) )
		} # Close INFO_PED loop
	} # Close "samp in KP_SAMPS" loop
	COMPILE <- list( PROBLEM_IPED,PROBLEM_PED, PROBLEM_IPED_RAND,PROBLEM_PED_RAND, HAVE_FAM_IPED,HAVE_FAM_PED, HAVE_FAM_IPED_RAND,HAVE_FAM_PED_RAND)
	names(COMPILE) <- c("IPED","PED","IPED_R","PED_R","T_IPED","T_PED","T_IPED_R","T_PED_R")
	return( COMPILE )
} # Close "REMOVE_RELATED" Function
PROBLEM_LIST.1 <- REMOVE_RELATED(KP_SAMPS)

length(unique(PROBLEM_LIST.1$IPED))
length(unique(PROBLEM_LIST.1$IPED_R))
length(unique(PROBLEM_LIST.1$PED))
length(unique(PROBLEM_LIST.1$PED_R))
length(intersect(PROBLEM_LIST.1$PED_R,PROBLEM_LIST.1$PED))
length(intersect(PROBLEM_LIST.1$IPED_R,PROBLEM_LIST.1$IPED))
length(union(PROBLEM_LIST.1$IPED,PROBLEM_LIST.1$PED))

## Try a few times
PROBLEM_LIST.2 <- REMOVE_RELATED(KP_SAMPS)
length(union(PROBLEM_LIST.2$IPED,PROBLEM_LIST.2$PED))
PROBLEM_LIST.3 <- REMOVE_RELATED(KP_SAMPS)
length(union(PROBLEM_LIST.3$IPED,PROBLEM_LIST.3$PED))
PROBLEM_LIST.4 <- REMOVE_RELATED(KP_SAMPS)
length(union(PROBLEM_LIST.4$IPED,PROBLEM_LIST.4$PED))
PROBLEM_LIST.5 <- REMOVE_RELATED(KP_SAMPS)
length(union(PROBLEM_LIST.5$IPED,PROBLEM_LIST.5$PED))

## Just for fun
length( Reduce( intersect, list(PROBLEM_LIST.5$PED,PROBLEM_LIST.4$PED) ))

## Remove Anyone that still has familial relationships w/ one another
RM_BC_FAM.1 <- union(PROBLEM_LIST.1$IPED,PROBLEM_LIST.1$PED)
KP_SAMPS.2.1 <- setdiff(KP_SAMPS, RM_BC_FAM.1 )
RM_BC_FAM.2 <- union(PROBLEM_LIST.2$IPED,PROBLEM_LIST.2$PED)
KP_SAMPS.2.2 <- setdiff(KP_SAMPS, RM_BC_FAM.2 )
RM_BC_FAM.3 <- union(PROBLEM_LIST.3$IPED,PROBLEM_LIST.3$PED)
KP_SAMPS.2.3 <- setdiff(KP_SAMPS, RM_BC_FAM.3 )
RM_BC_FAM.4 <- union(PROBLEM_LIST.4$IPED,PROBLEM_LIST.4$PED)
KP_SAMPS.2.4 <- setdiff(KP_SAMPS, RM_BC_FAM.4 )
RM_BC_FAM.5 <- union(PROBLEM_LIST.5$IPED,PROBLEM_LIST.5$PED)
KP_SAMPS.2.5 <- setdiff(KP_SAMPS, RM_BC_FAM.5 )

# ## Try checking for family again w/ new KP_SAMPS.2.X lists
# PROBLEM_LIST.2.1 <- REMOVE_RELATED(KP_SAMPS.2.1)
# length(union(PROBLEM_LIST.2.1$IPED,PROBLEM_LIST.2.1$PED))
# PROBLEM_LIST.2.2 <- REMOVE_RELATED(KP_SAMPS.2.2)
# length(union(PROBLEM_LIST.2.2$IPED,PROBLEM_LIST.2.2$PED))
# PROBLEM_LIST.2.3 <- REMOVE_RELATED(KP_SAMPS.2.3)
# length(union(PROBLEM_LIST.2.3$IPED,PROBLEM_LIST.2.3$PED))
# PROBLEM_LIST.2.4 <- REMOVE_RELATED(KP_SAMPS.2.4)
# length(union(PROBLEM_LIST.2.4$IPED,PROBLEM_LIST.2.4$PED))
# PROBLEM_LIST.2.5 <- REMOVE_RELATED(KP_SAMPS.2.5)
# length(union(PROBLEM_LIST.2.5$IPED,PROBLEM_LIST.2.5$PED))

######################################################
## GET SAMPLES BY ANCESTRY ###########################
######################################################

## Which Ancestry Options do I have?
 # http://www.1000genomes.org/faq/which-populations-are-part-your-study
 # SAN = SAS = GIH, PJL, BEB, STU, ITU
 # AMR = AMR = MXL, PUR, CLM, PEL
 # AFR = AFR = YRI, LWK, GWD, MSL, ESN, ASW, ACB
 # EUR = EUR = CEU, TSI, FIN, GBR, IBS
 # ASN = EAS = CHB, JPT, CHS, CDX, KHV
SAMPS_SAN <- as.character( KEY$sample[which(KEY$super_pop=="SAN")] )
SAMPS_ASN <- as.character( KEY$sample[which(KEY$super_pop=="ASN")] )
SAMPS_AMR <- as.character( KEY$sample[which(KEY$super_pop=="AMR")] )
SAMPS_AFR <- as.character( KEY$sample[which(KEY$super_pop=="AFR")] )
SAMPS_EUR <- as.character( KEY$sample[which(KEY$super_pop=="EUR")] )

######################################################
## GET FINAL 1KG SAMPLE LISTS/VCF ####################
######################################################

## Keep + Europeans
KP_INT_1KG.1 <- intersect( SAMPS_EUR, KP_SAMPS.2.1 )
KP_INT_1KG.2 <- intersect( SAMPS_EUR, KP_SAMPS.2.2 )
KP_INT_1KG.3 <- intersect( SAMPS_EUR, KP_SAMPS.2.3 )
KP_INT_1KG.4 <- intersect( SAMPS_EUR, KP_SAMPS.2.4 )
KP_INT_1KG.5 <- intersect( SAMPS_EUR, KP_SAMPS.2.5 )
length( Reduce( intersect, list(KP_INT_1KG.1,KP_INT_1KG.2,KP_INT_1KG.3,KP_INT_1KG.4,KP_INT_1KG.5) ))
 # Use KP_INT_1KG.5...2 extra people...most are the same b/n lists anyway

## Pull out Only Unrelated
X <- 1:12
VCF_EUR.C <- VCF_EUR[ ,c( 1:9,which(colnames(VCF_EUR) %in% KP_INT_1KG.5) ) ]
VCF_EUR.C[X,X]

############################################################################################################
## GET ALLELE FREQUENCIES FOR TBL AND 1KG COHORTS ##########################################################
############################################################################################################

######################################################
## FILTER/FIX GENOTYPES ##############################
######################################################

## Function to Switch Strand if there's a mismatch
STRAND <- function(GENOTYPE) {
	if (GENOTYPE=="A/A") { OUTPUT <- "T/T" }
	if (GENOTYPE=="A/C") { OUTPUT <- "T/G" }
	if (GENOTYPE=="A/G") { OUTPUT <- "T/C" }
	if (GENOTYPE=="A/T") { OUTPUT <- "T/A" }
	if (GENOTYPE=="C/A") { OUTPUT <- "G/T" }
	if (GENOTYPE=="C/C") { OUTPUT <- "G/G" }
	if (GENOTYPE=="C/G") { OUTPUT <- "G/C" }
	if (GENOTYPE=="C/T") { OUTPUT <- "G/A" }
	if (GENOTYPE=="G/A") { OUTPUT <- "C/T" }
	if (GENOTYPE=="G/C") { OUTPUT <- "C/G" }
	if (GENOTYPE=="G/G") { OUTPUT <- "C/C" }
	if (GENOTYPE=="G/T") { OUTPUT <- "C/A" }
	if (GENOTYPE=="T/A") { OUTPUT <- "A/T" }
	if (GENOTYPE=="T/C") { OUTPUT <- "A/G" }
	if (GENOTYPE=="T/G") { OUTPUT <- "A/C" }
	if (GENOTYPE=="T/T") { OUTPUT <- "A/A" }
	return(OUTPUT)
}

## Get list of Variants
VARS_1KG <- as.character( VCF_EUR.C$ID )
VARS_TBL <- colnames( GENO.C )[2:ncol(GENO.C)]
VARS_BOTH <- intersect( VARS_1KG, VARS_TBL )
VARS_1KG.only <- setdiff( VARS_1KG, VARS_TBL )
VARS_TBL.only <- setdiff( VARS_TBL, VARS_1KG )

## Filter for Variants present in Both data sets
 # TBL Data
WHICH.TBL <- which( colnames(GENO.C) %in% VARS_BOTH )
GENO.B <- GENO.C[, c(1, WHICH.TBL)]
 # 1KG Data
WHICH.1KG <- which( VCF_EUR.C$ID %in% VARS_BOTH )
VCF_EUR.B <- VCF_EUR.C[ which(VCF_EUR.C$ID %in% VARS_BOTH), ]

## Determine if TBL Alleles Match Ref from VCF
COMPILE_ALLELES <- array(, c(length(VARS_BOTH),3) )
rownames(COMPILE_ALLELES) <- VARS_BOTH ; colnames(COMPILE_ALLELES) <- c("TBL","1KG","MATCH")
COMPILE_ALLELES.F <- array(, c(length(VARS_BOTH),3) )
rownames(COMPILE_ALLELES.F) <- VARS_BOTH ; colnames(COMPILE_ALLELES.F) <- c("TBL","1KG","MATCH")
FLIPPED <- c()
for ( v in 1:length(VARS_BOTH) ) {
	var <- VARS_BOTH[v]
	WHICH_ROW <- which( VCF_EUR.B$ID==var )
	WHICH_COL <- which( colnames(GENO.B)==var )
	## Determine if Alleles/Strand Match(es)
	 # Get all alleles that show up in TBL data
	ALLELES_TBL.1 <- unique( c( sapply( strsplit(as.character(GENO.B[,WHICH_COL]),"/"), "[", 1:2 ), recursive=T) )
	ALLELES_TBL <- sort( ALLELES_TBL.1[which(ALLELES_TBL.1!="-")] )
	 # Get all alleles in 1KG data
	ALLELES_1KG <- sort( c( as.character(VCF_EUR.B$REF[WHICH_ROW]), as.character(VCF_EUR.B$ALT[WHICH_ROW]) ) )
	 # Find number of alleles that are on both
	ALLELES_MATCH_COUNT <- length( intersect( ALLELES_TBL, ALLELES_1KG ) )
	COMPILE_ALLELES[v,] <- c( paste(ALLELES_TBL,collapse="/"), paste(ALLELES_1KG,collapse="/"), ALLELES_MATCH_COUNT )
	## If Alleles/Strands don't match, flip it
	if ( ALLELES_MATCH_COUNT < 2 & length(ALLELES_TBL)==2 ) {
		NEW_ALLELES_TBL <- STRAND(paste(ALLELES_TBL,collapse="/"))
		COMPILE_ALLELES.F[v,] <- c( NEW_ALLELES_TBL, paste(ALLELES_1KG,collapse="/"), ALLELES_MATCH_COUNT )
		FLIPPED <- c( FLIPPED, var )
	}else{ COMPILE_ALLELES.F[v,] <- c( paste(ALLELES_TBL,collapse="/"), paste(ALLELES_1KG,collapse="/"), ALLELES_MATCH_COUNT )}	
}
COMPILE_ALLELES[FLIPPED,]

## Go through Rows in TBL Data that had Strand Issues
GENO.F <- GENO.B
for ( col in 1:ncol(GENO.F) ) { GENO.F[,col] <- as.character(GENO.F[,col]) }
HAD_STRAND_ISSUE <- which( COMPILE_ALLELES[,"MATCH"]==0 )
for ( v in HAD_STRAND_ISSUE ) {
	var <- rownames(COMPILE_ALLELES)[v]
	WHICH_COL <- which( colnames(GENO.F)==var )
	for ( row in 1:nrow(GENO.F) ) {
		GENO.F[row,WHICH_COL] <- STRAND( GENO.F[row,WHICH_COL] )
	}
}

######################################################
## 1KG ALLELE FREQUENCIES ############################
######################################################

## Go through each Person (column) and Change the Genotypes to Alt Allele Count
VCF_EUR.AAC <- VCF_EUR.B
for ( column in 10:ncol(VCF_EUR.AAC) ) {
	VCF_EUR.AAC[,column] <- as.character( VCF_EUR.AAC[,column] )
	VCF_EUR.AAC[which(VCF_EUR.AAC[,column]=="0|0"),column] <- 0
	VCF_EUR.AAC[which(VCF_EUR.AAC[,column]=="1|0"),column] <- 1
	VCF_EUR.AAC[which(VCF_EUR.AAC[,column]=="0|1"),column] <- 1
	VCF_EUR.AAC[which(VCF_EUR.AAC[,column]=="1|1"),column] <- 2
}
VCF_EUR.B[X,X]
VCF_EUR.AAC[X,X]

## Calculate Alternate Allele Frequencies for Population for each SNP
POP_1KG_AAC <- rowSums( data.matrix( VCF_EUR.AAC[,10:ncol(VCF_EUR.AAC)] ) )
POP_1KG_AAF <- POP_1KG_AAC / ( 2*(ncol(VCF_EUR.AAC)-9) )
names(POP_1KG_AAC) <- names(POP_1KG_AAF) <- VCF_EUR.AAC$ID

## Calculate Genotype Frequencies for Population for each SNP
POP_1KG_GC <- t( sapply( apply( VCF_EUR.AAC[,10:ncol(VCF_EUR.AAC)], 1, table ), "[", 1:3 ) )
POP_1KG_GC[which(is.na(POP_1KG_GC))] <- 0
POP_1KG_GF <- POP_1KG_GC / ( (ncol(VCF_EUR.AAC)-9) )
rownames(POP_1KG_GC) <- rownames(POP_1KG_GF) <- VCF_EUR.AAC$ID

######################################################
## TBL ALLELE FREQUENCIES ############################
######################################################

## Go through each SNP (column) and Change the Genotypes to Alt Allele Count
GENO.AAC <- GENO.F
for ( column in 2:ncol(GENO.AAC) ) {
	var <- colnames(GENO.AAC)[column]
	REF_ALL <- as.character( VCF_EUR.AAC$REF[which(VCF_EUR.AAC$ID==var)] )
	ALT_ALL <- as.character( VCF_EUR.AAC$ALT[which(VCF_EUR.AAC$ID==var)] )
	COL_SPLIT <- t(sapply( strsplit( as.character(GENO.AAC[,column]), "/" ), "[", 1:2 ))
	ALT_ALL_COUNT <- apply( COL_SPLIT, 1, function(X) length(which(X==ALT_ALL)) )
	GENO.AAC[,column] <- ALT_ALL_COUNT
}

## Calculate Alternate Allele Frequencies for Population for each SNP
POP_TBL_AAC <- colSums( data.matrix( GENO.AAC[,2:ncol(GENO.AAC)] ) )
POP_TBL_AAF <- POP_TBL_AAC / ( 2*(nrow(GENO.AAC)-9) )
names(POP_TBL_AAC) <- names(POP_TBL_AAF) <- colnames(GENO.AAC)[2:ncol(GENO.AAC)]

## Calculate Genotype Frequencies for Population for each SNP
POP_TBL_GC.LIST <- apply( GENO.AAC[,2:ncol(GENO.AAC)], 2, table )
POP_TBL_GC <- array(0, c(length(POP_TBL_GC.LIST),3) ) ; colnames(POP_TBL_GC) <- 0:2
for ( var in 1:length(POP_TBL_GC.LIST) ) { POP_TBL_GC[var,names(POP_TBL_GC.LIST[[var]])] <- POP_TBL_GC.LIST[[var]] }
POP_TBL_GF <- POP_TBL_GC / ( (nrow(GENO.AAC)) )
rownames(POP_TBL_GC) <- rownames(POP_TBL_GF) <- colnames(GENO.AAC)[2:ncol(GENO.AAC)]

######################################################
## COMBINE FREQUENCY TABLES ##########################
######################################################

## Calculate Alternate Allele Frequencies for Population for each SNP
POP_1KG <- data.frame( SNP=names(POP_1KG_AAC), AC=POP_1KG_AAC, AF=POP_1KG_AAF, GC=POP_1KG_GC, GF=POP_1KG_GF )
colnames(POP_1KG)[2:ncol(POP_1KG)] <- paste(colnames(POP_1KG)[2:ncol(POP_1KG)],"1KG",sep="_")
POP_TBL <- data.frame( SNP=names(POP_TBL_AAC), AC=POP_TBL_AAC, AF=POP_TBL_AAF, GC=POP_TBL_GC, GF=POP_TBL_GF )
colnames(POP_TBL)[2:ncol(POP_TBL)] <- paste(colnames(POP_TBL)[2:ncol(POP_TBL)],"TBL",sep="_")

## Merge Tables (hopefully by rowname)
POP_BOTH <- merge( POP_1KG, POP_TBL, by="SNP" )
VCF_ALLELE_INFO <- VCF_EUR.AAC[,c("ID","REF","ALT")]
POP <- merge( VCF_ALLELE_INFO, POP_BOTH, by.x="ID",by.y="SNP")
rownames(POP) <- POP$ID

## Which Variants are G<->C or A<->T (could be strand issues)
TROUBLE <- Reduce( union, list(which(POP$REF=="C" & POP$ALT=="G"),which(POP$REF=="G" & POP$ALT=="C"),which(POP$REF=="T" & POP$ALT=="A"),which(POP$REF=="A" & POP$ALT=="T")) )
TROUBLE_SNPS <- POP$ID[ TROUBLE ]

## Plot correlation b/n Allele Frequency
png( paste(PathToOut,"Scatter_Allele_Freq.png",sep="/"), height=1200,width=1200, pointsize=30)
plot( 0,0, type="n", xlab="1KG Allele Freq",ylab="TBL Allele Freq", main="MAF - TBL vs 1KG",xlim=c(0,1),ylim=c(0,1))
abline( h=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
abline( v=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
abline( 0,1, lty=2, col="black", lwd=2 )
COLS <- rep("dodgerblue2",nrow(POP)) ; COLS[TROUBLE] <- "firebrick2"
points(POP$AF_1KG,POP$AF_TBL, col=COLS, pch=20 )
legend( 0,1, col=c("dodgerblue2","firebrick2"), pch=20, legend=c("Confirmed","Potential Strand Issue") )
dev.off()

######################################################
## RUN TEST ON AC & GC ###############################
######################################################

## Sample Sizes for TBL/1KG
SIZE_TBL <- median(rowSums(POP_TBL_GC))
SIZE_1KG <- median(rowSums(POP_1KG_GC))

## Test all Variants for Allele/Genotype Counts
p_GC <- p_AC <- numeric(length(VARS_BOTH))
names(p_GC) <- names(p_AC) <- VARS_BOTH
for ( v in 1:length(VARS_BOTH) ) {
	var <- VARS_BOTH[v]# as.character(POP$ID[v])
	## Genotype Count
	GC_1KG <- as.numeric( POP[var,c("GC.0_1KG","GC.1_1KG","GC.2_1KG")] )
	GC_TBL <- as.numeric( POP[var,c("GC.0_TBL","GC.1_TBL","GC.2_TBL")] )
	GC_COMB <- data.matrix( rbind( GC_1KG, GC_TBL ) )
	p_GC[v] <- fisher.test(GC_COMB)$p.value
	## Allele Count
	AC_1KG <- as.numeric( POP[var,"AC_1KG"] )
	AC_TBL <- as.numeric( POP[var,"AC_TBL"] )
	AC_COMB <- rbind( c(AC_1KG,2*SIZE_1KG-AC_1KG), c(AC_TBL,2*SIZE_TBL-AC_TBL) )
	p_AC[v] <- fisher.test(AC_COMB)$p.value
}
P_VALS <- data.frame( GC=p_GC, AC=p_AC, TROUB=rep(1,length(p_GC)) )
P_VALS$TROUB[ which(rownames(P_VALS) %in% POP$ID[TROUBLE]) ] <- 2
head(P_VALS)
P_VALS[which(P_VALS$TROUB==2),]

## Plot P-Values
png( paste(PathToOut,"SNP_Pvals.png",sep="/"), height=1200,width=2600, pointsize=26)
plot( 0,0,type="n", xlim=c(1,nrow(P_VALS)), ylim=c(0,-log10(min(P_VALS))), xlab="", ylab="-log10(p)", main="TBL vs 1KG Frequencies", xaxt="n" )
axis( 1, at=1:nrow(P_VALS), labels=rownames(P_VALS), las=2 )
abline( h=seq(0,10,1), lty=2, lwd=1, col="grey50")
abline( h=-log10(.05/nrow(P_VALS)), lty=2, lwd=3, col="chartreuse2" )
legend( .6*nrow(P_VALS), -log10(min(P_VALS)), pch=c(15,16,18), legend=c("Genotype Count","Allele Count","Strand Issue?"), col=c(rep("dodgerblue2",2),"firebrick2") )
points( 1:nrow(P_VALS), -log10(P_VALS$GC), col=c("dodgerblue2","firebrick2")[P_VALS$TROUB], pch=15)
points( 1:nrow(P_VALS), -log10(P_VALS$AC), col=c("dodgerblue2","firebrick2")[P_VALS$TROUB], pch=16)
dev.off()

## Plot P-Values for AC vs GC
LIM <- c( 0, ceiling(-log10(min(P_VALS[,1:2]))) )
png( paste(PathToOut,"SNP_Pvals_AvG.png",sep="/"), height=1200,width=1200, pointsize=30)
plot( 0,0, type="n", xlim=LIM, ylim=LIM, main="P-Values - Allele vs Genotype Counts", xlab="Allele Count [-log10(p)]", ylab="Genotype Count [-log10(p)]")
abline( 0,1, lty=1, lwd=2, col="black" )
abline( h=seq(0,10,1), lty=2, lwd=1, col="grey50" )
abline( v=seq(0,10,1), lty=2, lwd=1, col="grey50" )
points( -log10(P_VALS$AC), -log10(P_VALS$GC), col=c("dodgerblue2","firebrick2")[P_VALS$TROUB], pch=20 )
dev.off()

######################################################
## IDENTIFY PHENOTYPES ASSOCIATED WITH VARS ##########
######################################################
 # PGR = Phenotype-Genotype Relationships

## Combine PGR info into one table
 # PGR_DIET, PGR_NUTR, PGR_EXER, PGR_META 
colnames(PGR_DIET) <- colnames(PGR_NUTR) <- colnames(PGR_EXER) <- colnames(PGR_META) <- c("PHENO","GENE","ID","RATE")
PGR <- rbind( data.frame(PGR_DIET,FUNC=rep("D",nrow(PGR_DIET))),
	data.frame(PGR_NUTR,FUNC=rep("N",nrow(PGR_NUTR))),
	data.frame(PGR_EXER,FUNC=rep("E",nrow(PGR_EXER))),
	data.frame(PGR_META,FUNC=rep("M",nrow(PGR_META))) )
PGR$ID <- gsub(" ","", as.character(PGR$ID) )

## Pull out SNPs w/ Decent P-Values
THRSH <- 1 # .001
WHICH_SNPS <- rownames(P_VALS)[ which( P_VALS$AC < THRSH | P_VALS$GC < THRSH ) ]
# WHICH_SNPS <- rownames(P_VALS)

## Merge Data into single table for Candidates
 # Include PGR Data
PGR_CANDS <- merge( x=data.frame(CAND=WHICH_SNPS), y=PGR, by.x="CAND",by.y="ID", all.x=T)
 # Include Population Frequency Data
CANDS.MG.1 <- merge( x=PGR_CANDS, y=POP, by.x="CAND", by.y="ID", all.x=T )
 # Include P-Values for Associaton Tests
CANDS.MG.2 <- merge( x=CANDS.MG.1, y=P_VALS, by.x="CAND", by.y="row.names", all.x=T )
 # Include Phenotype Pooling
CANDS.MG.3 <- merge( x=CANDS.MG.2, y=POOL, by.x="CAND",by.y="ID", all.x=T )
# CANDS.MG.3b <- array( ,c( 0,sum(ncol(CANDS.MG.2),ncol(POOL)) ) )
# colnames(CANDS.MG.3b) <- c( colnames(CANDS.MG.2), colnames(POOL) )
# for ( i in 1:nrow(CANDS.MG.2) ) {
# 	id <- as.character( CANDS.MG.2[i,"CAND"] )
# 	WHICH_ROW.pool <- grep( id, POOL$ID )
# 	if ( length(WHICH_ROW.pool)>0 ) {
# 		NEW_LINE <- c( CANDS.MG.2[i,], POOL[WHICH_ROW.pool,] )
# 	}else{ NEW_LINE <- c( CANDS.MG.2[i,], rep(NA,ncol(POOL)) ) }
# 	CANDS.MG.3b <- rbind( CANDS.MG.3b, NEW_LINE )
# }

CANDS.MG <- CANDS.MG.3

######################################################
## TEST POOLED GENOTYPE FREQUENCIES ##################
######################################################

BASE_FLIP <- function(base) {
	if ( base=="A" ) { out <- "T" }
	if ( base=="C" ) { out <- "G" }
	if ( base=="G" ) { out <- "C" }
	if ( base=="T" ) { out <- "A" }
	return(out)
}

## Loop through Variants
P.pool <- numeric( nrow(CANDS.MG) )
GT.pool <- list()
for ( i in 1:nrow(CANDS.MG) ) {
	id <- as.character( CANDS.MG[i,"CAND"] )
	GT_A <- as.character( CANDS.MG[i,"GT_A"] )
	GT_B <- as.character( CANDS.MG[i,"GT_B"] )
	TBL_GTS <- paste( GT_A, GT_B, sep="_" )
	if ( !grepl("NA",TBL_GTS) & !grepl("GT_B",TBL_GTS) ) {
		# Get Ref/Alt Alleles
		REF_ALLELE <- CANDS.MG[i,"REF"]
		ALT_ALLELE <- CANDS.MG[i,"ALT"]
		# Flip if necessary
		if ( !grepl(REF_ALLELE,TBL_GTS) & !grepl(ALT_ALLELE,TBL_GTS) ) {
			REF_ALLELE <- BASE_FLIP(REF_ALLELE)
			ALT_ALLELE <- BASE_FLIP(ALT_ALLELE)
		}
		# Add Genotype Frequencies for Pool
		GT_0 <- paste( REF_ALLELE, REF_ALLELE, sep="" )
		GT_1a <- paste( REF_ALLELE, ALT_ALLELE, sep="" )
		GT_1b <- paste( ALT_ALLELE, REF_ALLELE, sep="" )
		GT_2 <- paste( ALT_ALLELE, ALT_ALLELE, sep="" )
		# GTs <- c( GT_0, GT_1a, GT_1b, GT_2 )
		GT_A.spl <- strsplit( GT_A, "_" )[[1]]
		GT_B.spl <- strsplit( GT_B, "_" )[[1]]
		if ( GT_0 %in% GT_A.spl ) {
			print("GO")
			if ( GT_1a %in% GT_A.spl | GT_1b %in% GT_A.spl ) {
				GT_A.pool.tbl <- sum( CANDS.MG[i,c("GC.0_TBL","GC.1_TBL")] )
				GT_A.pool.1kg <- sum( CANDS.MG[i,c("GC.0_1KG","GC.1_1KG")] )
				GT_B.pool.tbl <- CANDS.MG[i,"GC.2_TBL"]
				GT_B.pool.1kg <- CANDS.MG[i,"GC.2_1KG"]
			}else{
				GT_A.pool.tbl <- CANDS.MG[i,"GC.0_TBL"]
				GT_A.pool.1kg <- CANDS.MG[i,"GC.0_1KG"]
				GT_B.pool.tbl <- sum( CANDS.MG[i,c("GC.1_TBL","GC.2_TBL")] )
				GT_B.pool.1kg <- sum( CANDS.MG[i,c("GC.1_1KG","GC.2_1KG")] )
			}
		}else{ # Means that GT_0 is NOT in GT_A.spl...so it must be in GT_B.spl
			print("OG")
			if ( GT_1a %in% GT_B.spl | GT_1b %in% GT_B.spl ) {
				GT_B.pool.tbl <- sum( CANDS.MG[i,c("GC.0_TBL","GC.1_TBL")] )
				GT_B.pool.1kg <- sum( CANDS.MG[i,c("GC.0_1KG","GC.1_1KG")] )
				GT_A.pool.tbl <- CANDS.MG[i,"GC.2_TBL"]
				GT_A.pool.1kg <- CANDS.MG[i,"GC.2_1KG"]
			}else{
				GT_B.pool.tbl <- CANDS.MG[i,"GC.0_TBL"]
				GT_B.pool.1kg <- CANDS.MG[i,"GC.0_1KG"]
				GT_A.pool.tbl <- sum( CANDS.MG[i,c("GC.1_TBL","GC.2_TBL")] )
				GT_A.pool.1kg <- sum( CANDS.MG[i,c("GC.1_1KG","GC.2_1KG")] )
			}
		}
		## Throw into Contingency Table
		GT.pool[[i]] <- array( c(GT_A.pool.tbl,GT_B.pool.tbl,GT_A.pool.1kg,GT_B.pool.1kg),c(2,2) )
		rownames(GT.pool[[i]]) <- c("GT_A","GT_B") ; colnames(GT.pool[[i]]) <- c("TBL","1KG")
		## Run Fisher Test
		P.pool[i] <- fisher.test( GT.pool[[i]] )$p.value
	}else{
		GT.pool[[i]] <- "No Pool Info"
		P.pool[i] <- 10
	}
	names(GT.pool)[i] <- id
	names(P.pool)[i] <- id
}

## Plots P-Values for various Tests
CANDS.UNQ <- data.frame( CANDS.MG, POOL=P.pool )
CANDS.UNQ <- CANDS.UNQ[ which(!duplicated(CANDS.UNQ[,"CAND"])), ]
COLS <- c("chartreuse2","dodgerblue2","slateblue2")
COLS.fac <- factor( CANDS.UNQ[,"TROUB"] )
XLIM <- c( 1,nrow(CANDS.UNQ) )
YLIM <- c( 0, -log10( min( CANDS.UNQ[,c("GC","AC","POOL")] ) ) )
png( paste(PathToOut,"SNP_Pvals_Pool.png",sep="/"), height=1200,width=2600, pointsize=26)
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, main="P-Values for Frequency Tests", ylab="-log10(p)", xaxt="n",xlab="")
axis(1, 1:nrow(CANDS.UNQ), label=CANDS.UNQ[,"CAND"], las=2 )
abline( h=seq( 0,100,1 ), lty=2, col="grey50" )
abline( h=-log10( .05/nrow(CANDS.UNQ) ), col="chocolate2", lty=2, lwd=3 )
abline( h=-log10( .05/nrow(CANDS.UNQ)/3 ), col="chocolate3", lty=2, lwd=3 )
legend( .6*nrow(CANDS.UNQ), YLIM[2], pch=c(20,4,3,18), legend=c("Genotype Count","Allele Count","Pooled Pheno","Strand Issue?"), col=c(COLS,"firebrick2"), pt.lwd=2, pt.cex=1.2 )
points( 1:nrow(CANDS.UNQ), -log10(CANDS.UNQ[,"GC"]), col=c(COLS[1],"firebrick2")[COLS.fac], pch=20, cex=1.2 )
points( 1:nrow(CANDS.UNQ), -log10(CANDS.UNQ[,"AC"]), col=c(COLS[2],"firebrick2")[COLS.fac], pch=4, lwd=2, cex=1.2 )
points( 1:nrow(CANDS.UNQ), -log10(CANDS.UNQ[,"POOL"]), col=c(COLS[3],"firebrick2")[COLS.fac], pch=3, lwd=2, cex=1.2 )
dev.off()

######################################################
## COMBINE INDIVIDUAL VARIANT FILES ##################
######################################################

## Get tables into merge-able format
 # TBL: GENO.AAC
AAC.TBL <- t( GENO.AAC[2:ncol(GENO.AAC)] )
colnames(AAC.TBL) <- GENO.AAC[,1]
 # 1KG: VCF_EUR.AAC
AAC.1KG <- VCF_EUR.AAC[ , 10:ncol(VCF_EUR.AAC) ]
rownames(AAC.1KG) <- VCF_EUR.AAC[,"ID"]

## Merge Tables
AAC.MRG.1 <- merge( x=AAC.TBL, y=AAC.1KG, by="row.names" )

## Get rid of "Trouble" variants
which_trouble <- which( AAC.MRG.1[,1] %in% TROUBLE_SNPS )
AAC.MRG.2 <- AAC.MRG.1[ -which_trouble, ]
AAC.MRG <- AAC.MRG.2[, 2:ncol(AAC.MRG.2) ]
rownames(AAC.MRG) <- AAC.MRG.2[,1]

## Merge Sample Names and Add Phenotype (i.e., Cohort Tag)
SAMPLE_NAMES.TBL <- as.character( GENO.AAC[,1] )
SAMPLE_NAMES.1KG <- colnames( VCF_EUR.AAC[, 10:ncol(VCF_EUR.AAC) ] )
SAMPLE_NAMES.ALL <- colnames( AAC.MRG )
AAC.PHENO <- data.frame( IID=SAMPLE_NAMES.ALL, FID=SAMPLE_NAMES.ALL, Pheno=rep(1,length(SAMPLE_NAMES.ALL)) )
AAC.PHENO$Pheno[ which(AAC.PHENO$IID %in% SAMPLE_NAMES.1KG) ] <- 2

## Convert AAC files to VCF format
which_snps <- which( VCF_EUR$ID %in% rownames(AAC.MRG) )
VCF.MRG.1 <- VCF_EUR.AAC[ which_snps, 1:9 ]
VCF.MRG.2 <- merge( x=VCF.MRG.1, y=AAC.MRG, by.x="ID",by.y="row.names" )
which_0 <- which(VCF.MRG.2==0, arr.ind=T) ; which_0 <- which_0[ which(which_0[,2]>=10), ]
which_1 <- which(VCF.MRG.2==1, arr.ind=T) ; which_1 <- which_1[ which(which_1[,2]>=10), ]
which_2 <- which(VCF.MRG.2==2, arr.ind=T) ; which_2 <- which_2[ which(which_2[,2]>=10), ]
VCF.MRG.2[which_0] <- "0/0"
VCF.MRG.2[which_1] <- "0/1"
VCF.MRG.2[which_2] <- "1/1"

## Write these files to a directory
write.table( VCF.MRG.2, file=paste(PathToOut,"/Merged_AAC.vcf",sep=""), sep="\t",row.names=F, col.names=T,quote=F )
write.table( AAC.PHENO, file=paste(PathToOut,"/Merged_Pheno.txt",sep=""), sep="\t",row.names=F, col.names=T,quote=F )

######################################################
## PRINCIPAL COMPONENT ANALYSIS ######################
######################################################

## Calculate PCs from Variants
SDS <- apply( AAC.MRG, 1, sd )
AAC.MRG.2 <- AAC.MRG[ -which(SDS==0), ]
PC.1 <- prcomp( t( data.matrix(AAC.MRG.2) ), scale=T, center=T )
COLS <- c("springgreen2","slateblue4")
pairs( PC.1$x[,1:5], col=COLS[factor(AAC.PHENO$Pheno)], pch="+" )
identical( rownames(PC.1$x), as.character(AAC.PHENO$IID) )

######################################################
## LASSO/ELASTIC NET REGRESSION ######################
######################################################

## Prep for 10X cross-validation
N_SAMPS <- length(SAMPLE_NAMES.ALL )
SAMP_ORDER <- sample( SAMPLE_NAMES.ALL, N_SAMPS, replace=F )
N_ITER <- 5
N_TEST <- ceiling( N_SAMPS / N_ITER )
L.bins <- c("BEST","SE","MIN","MAX")

BETA <- LAMBDA <- PRED <- PRED.CL <- list()
## Loop through N Iterations
start_time <- proc.time()
for ( i in 1:N_ITER ) {
	iter <- paste( "I",i,sep="_" )
	print(paste("## Iteration",i,"of",N_ITER,"-",round(proc.time()-start_time,3)[3] ))
	## Sample test/training sets
	which_test <- (i-1)*N_TEST + 1:N_TEST
	which_test <- which_test[ which(which_test <= N_SAMPS) ]
	which_train <- setdiff( 1:N_SAMPS, which_test )
	test_set <- SAMP_ORDER[ which_test ]
	train_set <- SAMP_ORDER[ which_train ]
	## Sample Phenotype/Genotype Data
	Y.train <- data.matrix( AAC.PHENO[ which(AAC.PHENO$IID %in% train_set), "Pheno" ] )
	Y.test <- data.matrix( AAC.PHENO[ which(AAC.PHENO$IID %in% test_set), "Pheno" ] )
	X.train <- t(data.matrix( AAC.MRG[, which(AAC.PHENO$IID %in% train_set) ] ))
	X.test <- t(data.matrix( AAC.MRG[, which(AAC.PHENO$IID %in% test_set) ] ))
	## Run Lasso Regression
	 # Set Fit Parameters
	INT <- F
	MEAS <- "class"
	 # Fit this Shiz
	print("Fitting Models")
	fit0 <- glmnet(x=X.train, y=Y.train, family="binomial", alpha=0, intercept=INT )
	fit.5 <- glmnet(x=X.train, y=Y.train, family="binomial", alpha=.5, intercept=INT )
	fit1 <- glmnet(x=X.train, y=Y.train, family="binomial", alpha=1, intercept=INT )
	fit.cv0 <- cv.glmnet(x=X.train, y=Y.train, family="binomial", alpha=0, intercept=INT, nfolds=20, type.measure=MEAS )
	fit.cv.5 <- cv.glmnet(x=X.train, y=Y.train, family="binomial", alpha=.5, intercept=INT, nfolds=20, type.measure=MEAS )
	fit.cv1 <- cv.glmnet(x=X.train, y=Y.train, family="binomial", alpha=1, intercept=INT, nfolds=20, type.measure=MEAS )
	## Pull Lambda Values
	L0 <- c( fit.cv0$lambda.min, fit.cv0$lambda.1se, min(fit.cv0$lambda), max(fit.cv0$lambda) )
	L.5 <- c( fit.cv.5$lambda.min, fit.cv.5$lambda.1se, min(fit.cv.5$lambda), max(fit.cv.5$lambda) )
	L1 <- c( fit.cv1$lambda.min, fit.cv1$lambda.1se, min(fit.cv1$lambda), max(fit.cv1$lambda) )
	L <- rbind( L0, L.5, L1 )
	colnames(L) <- c("Best","SE","Min","Max")
	rownames(L) <- paste("Alpha=",c(0,.5,1),sep="")
	 # Plot Fits w/ Different Alpha Values
	COLS <- c("springgreen2","steelblue2","slateblue2")
	YLIM0 <- c( min(fit.cv0$cvlo), max(fit.cv0$cvup) )
	YLIM.5 <- c( min(fit.cv.5$cvlo), max(fit.cv.5$cvup) )
	YLIM1 <- c( min(fit.cv1$cvlo), max(fit.cv1$cvup) )
	L.COLS0 <- paste("springgreen",1:4,sep="") # c("firebrick2","gold2","chartreuse2","dodgerblue2")
	L.COLS.5 <- paste("steelblue",1:4,sep="")
	L.COLS1 <- paste("slateblue",1:4,sep="")
	# XLIM <- log( range( c(fit.cv0$lambda, fit.cv.5$lambda, fit.cv1$lambda) ) )
	# YLIM <- extendrange( c(fit.cv0$cvm,fit.cv.5$cvm,fit.cv1$cvm) )
	png( paste(PathToOut,"/CVfit_",i,"_AB.png",sep=""), height=1600,width=2400, pointsize=32)
	par( mfrow=c(2,3) )
	par( mar=c(5,5,5,3) )
	# plot( fit.cv0, main="Alpha=0 (Ridge)" )
	# plot( fit.cv.5, main="Alpha=.5" )
	# plot( fit.cv1, main="Alpha=1 (LASSO)" )
	plot( log(fit.cv0$lambda),fit.cv0$cvm, main="Alpha=0 (Ridge)", pch=20,col=COLS[1],xlab="log(Lambda)",ylab=fit.cv0$name, ylim=YLIM0)
	arrows( log(fit.cv0$lambda),fit.cv0$cvup,log(fit.cv0$lambda),fit.cv0$cvlo, code=3,angle=90,length=.05 )
	abline( h=seq(0,1,.1), lty=c(rep(2,5),1,rep(2,5)), col=c(rep("grey50",5),"grey20",rep("grey50",5)))
	abline( v=log(L0), col=L.COLS0, lwd=2 )
	plot( log(fit.cv.5$lambda),fit.cv.5$cvm, main="Alpha=.5", pch=20,col=COLS[2],xlab="log(Lambda)",ylab=fit.cv.5$name, ylim=YLIM.5)
	arrows( log(fit.cv.5$lambda),fit.cv.5$cvup,log(fit.cv.5$lambda),fit.cv.5$cvlo, code=3,angle=90,length=.05 )
	abline( h=seq(0,1,.1), lty=c(rep(2,5),1,rep(2,5)), col=c(rep("grey50",5),"grey20",rep("grey50",5)))
	abline( v=log(L.5), col=L.COLS.5, lwd=2 )
	plot( log(fit.cv1$lambda),fit.cv1$cvm, main="Alpha=1 (LASSO)", pch=20,col=COLS[3],xlab="log(Lambda)",ylab=fit.cv1$name, ylim=YLIM1)
	arrows( log(fit.cv1$lambda),fit.cv1$cvup,log(fit.cv1$lambda),fit.cv1$cvlo, code=3,angle=90,length=.05 )
	abline( h=seq(0,1,.1), lty=c(rep(2,5),1,rep(2,5)), col=c(rep("grey50",5),"grey20",rep("grey50",5)))
	abline( v=log(L1), col=L.COLS1, lwd=2 )
	 # Plot Fits
	plot( fit0, xvar="lambda", main="Alpha=0 (Ridge)" ) ; abline( v=log(L0), col=L.COLS0, lwd=2 )
	plot( fit.5, xvar="lambda", main="Alpha=.5" ) ; abline( v=log(L.5), col=L.COLS.5, lwd=2 )
	plot( fit1, xvar="lambda", main="Alpha=1 (LASSO)" ) ; abline( v=log(L1), col=L.COLS1, lwd=2 )
	dev.off()
	 # Barplot of Lambda Values
	YLIM <- range(log(L))+c(0,4)
	png( paste(PathToOut,"/CVfit_",i,"_Lambdas.png",sep=""), height=1000,width=1200, pointsize=32)
	barplot( log(L), col=COLS, beside=T, legend=T, ylab="log(Lambda)", main="Lambda Values for Fits", ylim=YLIM )
	abline( h=seq(floor(YLIM[1]),ceiling(YLIM[2]),1), lty=2, col="grey50" )
	barplot( log(L), col=COLS, beside=T, legend=T, ylab="log(Lambda)", main="Lambda Values for Fits", ylim=YLIM, add=T )
	dev.off()
	## Pull Beta Values
	B0 <- data.matrix( coef(fit.cv0, s=L0[1] ) )
	B.5 <- data.matrix( coef(fit.cv.5, s=L.5[1] ) )
	B1 <- data.matrix( coef(fit.cv1, s=L1[1] ) )
	B <- data.frame( B0, B.5, B1 ) ; colnames(B) <- c("A0","A.5","A1")
	## Calculate Predicted Values
	 # Class Probability
	pred0 <- data.matrix( predict(fit.cv0, newx=X.test, s=L0[1], type="response" ) )
	pred.5 <- data.matrix( predict(fit.cv.5, newx=X.test, s=L.5[1], type="response" ) )
	pred1 <- data.matrix( predict(fit.cv1, newx=X.test, s=L1[1], type="response" ) )
	pred <- data.frame( OBS=Y.test, A0=pred0, A.5=pred.5, A1=pred1 )
	colnames(pred) <- c("OBS","A0","A.5","A1")
	COLS.4 <- gsub("2","4",COLS)
	png( paste(PathToOut,"/CVfit_",i,"_Class_Prob.png",sep=""), height=1000,width=1600, pointsize=32)
	plot( 0,0,type="n", xlim=c(0.5,3.5),ylim=c(0,1), main="Class Probability Distributions",xlab="Alpha Value",ylab="Class Probability", xaxt="n" )
	abline( h=seq(0,1,.1), lty=2, col="grey50" )
	boxplot( pred[,2:4], col=COLS, main="Class Probability Distributions",xlab="Alpha Value",ylab="Class Probability", pch="", add=T )
	for ( col in 2:4 ) { points( jitter(rep(col-1,nrow(pred))),pred[,col], col=COLS.4[col-1], pch=c("o","+")[factor(pred[,1])] ) }
	dev.off()
	# BRKS <- seq(0,1,.05)
	# plot( 0,0,type="n", xlim=c(0,1), ylim=c(0,20), xlab="Class Probability",ylab="Frequency" )
	# hist( pred[,2], breaks=BRKS, density=35, col=COLS[1], angle=30, add=T )
	# hist( pred[,3], breaks=BRKS, density=35, col=COLS[2], angle=-30, add=T )
	# hist( pred[,4], breaks=BRKS, density=35, col=COLS[3], angle=60, add=T )
	 # Class
	pred0.cl <- data.matrix( predict(fit.cv0, newx=X.test, s=L0[1], type="class" ) )
	pred.5.cl <- data.matrix( predict(fit.cv.5, newx=X.test, s=L.5[1], type="class" ) )
	pred1.cl <- data.matrix( predict(fit.cv1, newx=X.test, s=L1[1], type="class" ) )
	pred.cl <- data.frame( OBS=Y.test, A0=pred0.cl, A.5=pred.5.cl, A1=pred1.cl )
	# pairs( pred )
			## Compile Outputs
	print("Compiling Data")
	LAMBDA[[iter]] <- L
	BETA[[iter]] <- B
	PRED[[iter]] <- pred
	PRED.CL[[iter]] <- pred.cl
}
# for (j in 1:N_ITER) { for ( i in 2:4) { print( table( PRED.CL[[j]][,c(1,i)] ) ) } }

## Compile Predicted Values
PRED.2 <- PRED.CL.2 <- array( ,c(0,4) )
for ( i in 1:N_ITER ) {
	PRED.2 <- rbind( PRED.2, PRED[[i]] )
	PRED.CL.2 <- rbind( PRED.CL.2, PRED.CL[[i]] )
}
colnames(PRED.2) <- colnames(PRED.CL.2) <- c("OBS","A0","A.5","A1")
for ( i in 2:4 ) { print( table( PRED.CL.2[,c(1,i)] ) ) }

## Plot Histograms of Class Probabilities
BRKS <- seq(0,1,.05)
COLS <- c("black","springgreen2","steelblue2","slateblue2")
par(mfrow=c(3,1))
for ( i in 2:4 ) { hist( PRED.2[,i], breaks=BRKS, col=COLS[i] ) }

######################################################
## POLY TRAIN (REG) ##################################
######################################################

## Prep for 10X cross-validation
N_ITER <- 5
N_TEST <- ceiling( N_SAMPS / N_ITER )
MOD_SIZE <- 30

CV_COMP <- list()
## Loop through N Iterations
for ( n in 1:N_ITER ) {
	iter <- paste( "I",n,sep="_" )
	## Sample test/training sets
	which_test <- (n-1)*N_TEST + 1:N_TEST
	which_test <- which_test[ which(which_test <= N_SAMPS) ]
	which_train <- setdiff( 1:N_SAMPS, which_test )
	test_set <- SAMP_ORDER[ which_test ]
	train_set <- SAMP_ORDER[ which_train ]
	## Sample Phenotype/Genotype Data
	Y.train <- data.matrix( AAC.PHENO[ which(AAC.PHENO$IID %in% train_set), "Pheno" ] )
	Y.test <- data.matrix( AAC.PHENO[ which(AAC.PHENO$IID %in% test_set), "Pheno" ] )
	X.train <- t(data.matrix( AAC.MRG[, which(AAC.PHENO$IID %in% train_set) ] ))
	X.test <- t(data.matrix( AAC.MRG[, which(AAC.PHENO$IID %in% test_set) ] ))

	## Loop through Predictors and add each one iteratively to model
	 # Then assess fit of model
	# TO_TEST.gt <- length(TO_TEST.cov) + 1:length(OR.gt.9)
	TO_TEST <- 1:ncol(X.train)
	TO_INCL <- c()
	TO_OMIT <- c()
	N_PREDS <- ncol(X.train)
	MOD_SIZE <- min( MOD_SIZE, N_PREDS )
	AUC.test <- list()
	AUC.test$TS <- AUC.test$TR <- array( , c(N_PREDS,MOD_SIZE) )
	rownames(AUC.test$TS) <- rownames(AUC.test$TR) <- colnames(X.train)[1:N_PREDS]
	AUC <- array( , dim=c(MOD_SIZE,2) )
	colnames(AUC) <- c("TR","TS")
	# for ( r in 1:nrow(COEF.all) ) {
	print(paste( "## Building Models -",n,"of",N_ITER,"##" ))
	start_time <- proc.time()
	for ( i in 1:MOD_SIZE ) {
	# for ( i in 1:5 ) {
		## Loop through all Predictors and Find Best on to ADD to model
		for ( r in TO_TEST ) {
			This_Pred <- colnames(X.train)[r]
			Which_Pred <- c( TO_INCL, This_Pred )
			## Specify Phenotype Values
			PHENO.TR.9 <- factor(Y.train[,1])
			PHENO.TS.9 <- factor(Y.test[,1])
			## Build Temporary Data Frame with These Covariates
			TEMP.TR.9 <- data.matrix( X.train[, Which_Pred ] ) ; colnames(TEMP.TR.9) <- Which_Pred
			TEMP.TS.9 <- data.matrix( X.test[, Which_Pred ] ) ; colnames(TEMP.TS.9) <- paste("TEMP.TR.9",Which_Pred,sep="")
			## Fit Model Using lm()
			MOD.TR.9 <- glm( PHENO.TR.9 ~ TEMP.TR.9, family=binomial(link=logit) )
			COEF.mod <- coef(MOD.TR.9)
			PRED.TR.9 <- COEF.mod[1] + TEMP.TR.9 %*% matrix(COEF.mod[-1])
			PRED.TS.9 <- COEF.mod[1] + TEMP.TS.9 %*% matrix(COEF.mod[-1])
			## Calculate Probability from Link Function
			 # p <- exp( b0 + b1X1 + ... + bnXn ) / ( 1 + exp( b0 + b1X1 + ... + bnXn ) )
			PROB.TR.9 <- exp(PRED.TR.9) / ( 1 + exp(PRED.TR.9) )
			PROB.TS.9 <- exp(PRED.TS.9) / ( 1 + exp(PRED.TS.9) )
			## Loop through different thresholds for calling Case/Control
			PROB.UNION <- sort( union( PROB.TR.9, PROB.TS.9 ) )
			PROB.RANGE <- range( PROB.UNION )
			PROB.N <- length(PROB.UNION)
			if ( PROB.N < 100 ) {
				THRESHOLDS <- ( c(0,PROB.UNION)+c(PROB.UNION,1) ) / 2
			}else{
				THRESHOLDS.1 <- ( c(0,PROB.UNION)+c(PROB.UNION,1) ) / 2
				which_thresh <- c( 1, floor(seq(2,PROB.N,length.out=100)), PROB.N )
				THRESHOLDS <- THRESHOLDS.1[which_thresh]
				# BIN_SIZE <- .005
				# THRESHOLDS.1 <- seq( PROB.RANGE[1]-BIN_SIZE,PROB.RANGE[2]+BIN_SIZE,.005 )
				# THRESHOLDS.2 <- seq( PROB.RANGE[1]-BIN_SIZE,PROB.RANGE[2]+BIN_SIZE, length.out=20 )
				# if (length(THRESHOLDS.1)>length(THRESHOLDS.2) ) {
				# 	THRESHOLDS <- THRESHOLDS.1
				# }else{ THRESHOLDS <- THRESHOLDS.2 }		
			}
			# THRESHOLDS <- seq( 0,1,.002 )
			N_THRSH <- length(THRESHOLDS)
			TPR.TR.9 <- TPR.TS.9 <- FPR.TR.9 <- FPR.TS.9 <- numeric(N_THRSH)
			for ( t in 1:N_THRSH ) {
				thresh <- THRESHOLDS[t]
				## Predict Cases/Controls
				BIN.TR.9 <- PROB.TR.9 > thresh
				BIN.TS.9 <- PROB.TS.9 > thresh
				## Calculate False Positive and True Postive Rates
				TPR.TR.9[t] <- length(which( PHENO.TR.9==2 & BIN.TR.9==T )) / length(which( PHENO.TR.9==2 ))
				TPR.TS.9[t] <- length(which( PHENO.TS.9==2 & BIN.TS.9==T )) / length(which( PHENO.TS.9==2 ))
				FPR.TR.9[t] <- length(which( PHENO.TR.9==1 & BIN.TR.9==T )) / length(which( PHENO.TR.9==1 ))
				FPR.TS.9[t] <- length(which( PHENO.TS.9==1 & BIN.TS.9==T )) / length(which( PHENO.TS.9==1 ))
			}
			if ( runif(1,0,1) < .001 ) {
				plot( 0,0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="FPR",ylab="TPR", main="ROC" )
				abline( 0,1, lty=1, col="black", lwd=2 )
				abline( h=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
				abline( v=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
				points( FPR.TR.9, TPR.TR.9, type="o", pch=20, col="chartreuse2", lwd=2 )
				points( FPR.TS.9, TPR.TS.9, type="o", pch=20, col="dodgerblue2", lwd=2 )
			}
			## Calculate Area Under Curve
			AUC.test$TR[r,i] <- sum( -( FPR.TR.9[2:N_THRSH]-FPR.TR.9[2:N_THRSH-1] ) * ( TPR.TR.9[2:N_THRSH]+TPR.TR.9[2:N_THRSH-1] ) /2 )
			AUC.test$TS[r,i] <- sum( -( FPR.TS.9[2:N_THRSH]-FPR.TS.9[2:N_THRSH-1] ) * ( TPR.TS.9[2:N_THRSH]+TPR.TS.9[2:N_THRSH-1] ) /2 )
		} # Close "r" Loop
		## Decide Which Model is the Best
		WHICH_TO_INC <- which.max( AUC.test$TR[,i] )
		Pred_To_Inc <- colnames(X.train)[WHICH_TO_INC]
		print(paste( "Include",WHICH_TO_INC,"-",Pred_To_Inc ))
		## Compile all the results for Plotting
		AUC[i,"TR"] <- AUC.test$TR[WHICH_TO_INC,i]
		AUC[i,"TS"] <- AUC.test$TS[WHICH_TO_INC,i]
		## Include, but don't Test, that Predictor from now on
		TO_INCL <- c( TO_INCL, Pred_To_Inc )
		TO_TEST <- setdiff( TO_TEST, WHICH_TO_INC )
		## Output Status
		if (i%%5==0) { 
			# print( tail(Which_Pred) )
			print(paste( "Finished",i,"of",MOD_SIZE,"-",round( (proc.time()-start_time)[3], 1) ))
		}
	} # Close "i" Loop
	rownames(AUC) <- TO_INCL
	colnames(AUC.test$TS) <- colnames(AUC.test$TR) <- TO_INCL
	CV_COMP[[iter]] <- AUC
	## Plot Latest Iteration Performance
	XLIM <- c(1,nrow(AUC))
	YLIM <- c(0,1)
	png( paste(PathToOut,"/Poly_Train_",n,"_9_AUC.png",sep=""), height=1000,width=1200, pointsize=24)
	plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, main="AUC vs Model", ylab="AUC", xlab="",xaxt="n")
	axis( 1, at=1:nrow(AUC), label=rownames(AUC), las=2)
	abline( h=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
	abline( h=c(0,.5,1), lty=1, col="grey20", lwd=1 )
	abline( v=seq(0,XLIM[2],5), lty=2, col="grey50", lwd=1 )
	points( 1:nrow(AUC), AUC[,1], pch=20,type="o", col="chartreuse2", lwd=3)
	points( 1:nrow(AUC), AUC[,2], pch=20,type="o", col="deepskyblue2", lwd=3)
	dev.off()
} # Close "n" Loop

for ( n in 1:N_ITER) {
	AUC.plot <- CV_COMP[[n]]
	XLIM <- c(1,nrow(AUC.plot))
	YLIM <- c(0,1)
	png( paste(PathToOut,"/Poly_Train_",n,"_9_AUC.png",sep=""), height=1200,width=1600, pointsize=24)
	plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, main="AUC vs Model", ylab="AUC", xlab="",xaxt="n")
	axis( 1, at=1:nrow(AUC.plot), label=rownames(AUC.plot), las=2)
	abline( h=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
	abline( h=c(0,.5,1), lty=1, col="grey20", lwd=1 )
	abline( v=seq(0,XLIM[2],5), lty=2, col="grey50", lwd=1 )
	points( 1:nrow(AUC.plot), AUC.plot[,1], pch=20,type="o", col="chartreuse2", lwd=3)
	points( 1:nrow(AUC.plot), AUC.plot[,2], pch=20,type="o", col="deepskyblue2", lwd=3)
	dev.off()
}

##################################################
## BINARY: TBL v 1KG
RUNS <- names(CV_COMP)

## Set up variables to compile
 # List of variants in each model
VARS.9 <- VARS.1 <- list()
 # TS/TR Improvement for inclusion of each variant
AUCUP.9 <- AUCUP.1 <- list()

## Loop through Models & Compile Stats
n <- 0
for ( iter in RUNS ) {
	n <- n + 1
	## Get variant IDs
	VARS.9[[iter]] <- rownames( CV_COMP[[iter]] )
	## Get change in AUC from inclusion of each variant
	len.9 <- length(VARS.9[[iter]])
	AUCUP.9[[iter]] <- array( , c(len.9,2) )
	colnames(AUCUP.9[[iter]]) <- c("TR","TS")
	rownames(AUCUP.9[[iter]]) <- VARS.9[[iter]]
	AUCUP.9[[iter]][,"TR"] <- CV_COMP[[iter]][ VARS.9[[iter]], "TR" ] - c(0.5,CV_COMP[[iter]][ 1:(len.9-1), "TR" ] )
	AUCUP.9[[iter]][,"TS"] <- CV_COMP[[iter]][ VARS.9[[iter]], "TS" ] - c(0.5,CV_COMP[[iter]][ 1:(len.9-1), "TS" ] )
}

ALL.VARS.9 <- Reduce( union, VARS.9 )
RANKS.9 <- WTS.9 <- array( , c(length(ALL.VARS.9),length(RUNS)) )
rownames(RANKS.9) <- rownames(WTS.9) <- ALL.VARS.9
colnames(RANKS.9) <- colnames(WTS.9) <- RUNS
for ( v in 1:length(ALL.VARS.9) ) {
	var <- ALL.VARS.9[v]
	for ( run in RUNS ) {
		if ( var %in% VARS.9[[run]] ) {
			RANKS.9[var,run] <- which( VARS.9[[run]] == var )
			WTS.9[var,run] <- AUCUP.9[[run]][var,"TS"]
		}
	}
}

 # Count how many models each variant appears in
COUNT.9 <- apply( RANKS.9, 1, function(x) length(which(!is.na(x))) )
 # Get weight based on R2 improvement
WT.9 <- apply( WTS.9, 1, sum, na.rm=T )
 # Get best rank across models for each variant
MIN.9 <- apply( RANKS.9, 1, min, na.rm=T )
 # Calculate the MEAN rank of each variant
MEAN.9 <- apply( RANKS.9, 1, mean, na.rm=T )
 # Take MEAN rank and divide by number of runs in which variant appears
 # This boosts variants that appear in multiple models by giving a lower rank STAT
STAT.9 <- MIN.9 / COUNT.9

## Compile into single table
COMP.9 <- data.frame( COUNT=COUNT.9, MEAN_RNK=round(MEAN.9,2), STAT=round(STAT.9,2), AUC_SUM=round(WT.9,2), RANKS.9 )
COMP.9 <- COMP.9[order(COMP.9[,"STAT"]),]
COMP.9.AUC <- COMP.9[order(COMP.9[,"AUC_SUM"],decreasing=T),]

## Pull out Top s SNPs
s <- 4
TOP_SNPS <- rownames(COMP.9)[1:s]
Y <- as.factor( AAC.PHENO[,"Pheno"] )
X <- t( data.matrix( AAC.MRG[TOP_SNPS,] ) )
XY <- data.frame( X, Y )
MOD.top <- glm( Y ~ ., data=XY, family=binomial(link=logit) )
COEF.top <- coef(MOD.top)
PRED.top <- COEF.top[1] + X %*% matrix(COEF.top[-1])
PROB.top <- exp(PRED.top) / ( 1 + exp(PRED.top) )
PROB.UNION <- sort( unique(PROB.top) )
PROB.RANGE <- range( PROB.UNION )
PROB.N <- length(PROB.UNION)
if ( PROB.N < 100 ) {
	THRESHOLDS <- ( c(0,PROB.UNION)+c(PROB.UNION,1) ) / 2
}else{
	THRESHOLDS.1 <- ( c(0,PROB.UNION)+c(PROB.UNION,1) ) / 2
	which_thresh <- c( 1, floor(seq(2,PROB.N,length.out=100)), PROB.N )
	THRESHOLDS <- THRESHOLDS.1[which_thresh]
}
N_THRSH <- length(THRESHOLDS)
TPR.top <- FPR.top <- numeric(N_THRSH)
for ( t in 1:N_THRSH ) {
	thresh <- THRESHOLDS[t]
	## Predict Cases/Controls
	BIN.top <- PROB.top > thresh
	## Calculate False Positive and True Postive Rates
	TPR.top[t] <- length(which( Y==2 & BIN.top==T )) / length(which( Y==2 ))
	FPR.top[t] <- length(which( Y==1 & BIN.top==T )) / length(which( Y==1 ))
}
plot( 0,0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="FPR",ylab="TPR", main="ROC" )
abline( 0,1, lty=1, col="black", lwd=2 )
abline( h=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
abline( v=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
points( FPR.top, TPR.top, type="o", pch=20, col="chartreuse2", lwd=2 )

######################################################
## POLY TRAIN (SCORE) ################################
######################################################

## Compile P-Values and Beta Values for Covariates and Genotypes
ORD.cov <- order( P.cov ) ; ORD.gt <- order( P.gt.1 )
P.all <- c( P.cov[ORD.cov], P.gt.1[ORD.gt] ) ; names(P.all) <- c( names(P.cov)[ORD.cov], names(P.gt.1)[ORD.gt] )
OR.all <- c( OR.cov[ORD.cov], OR.gt.1[ORD.gt] ) ; names(OR.all) <- c( names(OR.cov)[ORD.cov], names(OR.gt.1)[ORD.gt] )
COEF.all <- data.frame( P.all, OR.all ) ; rownames(COEF.all) <- names(OR.all)
COEF.int <- COEF.all[ which(rownames(COEF.all)=="(Intercept)"), ]
COEF.all <- COEF.all[ -which(rownames(COEF.all)=="(Intercept)"), ]
COEF.all$OR.all[(length(OR.cov)):nrow(COEF.all)] <- c(-1,1)[factor(COEF.all$OR.all[(length(OR.cov)):nrow(COEF.all)]>0)]

## Prep for 10X cross-validation
N_ITER <- 5
N_TEST <- ceiling( N_SAMPS / N_ITER )
MOD_SIZE <- 30

CV_COMP <- list()
## Loop through N Iterations
for ( n in 1:N_ITER ) {
	iter <- paste( "I",n,sep="_" )
	## Sample test/training sets
	which_test <- (n-1)*N_TEST + 1:N_TEST
	which_test <- which_test[ which(which_test <= N_SAMPS) ]
	which_train <- setdiff( 1:N_SAMPS, which_test )
	test_set <- SAMP_ORDER[ which_test ]
	train_set <- SAMP_ORDER[ which_train ]
	## Sample Phenotype/Genotype Data
	Y.train <- data.matrix( AAC.PHENO[ which(AAC.PHENO$IID %in% train_set), "Pheno" ] )
	Y.test <- data.matrix( AAC.PHENO[ which(AAC.PHENO$IID %in% test_set), "Pheno" ] )
	X.train <- t(data.matrix( AAC.MRG[, which(AAC.PHENO$IID %in% train_set) ] ))
	X.test <- t(data.matrix( AAC.MRG[, which(AAC.PHENO$IID %in% test_set) ] ))

	## Loop through Predictors and add each one iteratively to model
	 # Then assess fit of model
	# TO_TEST.gt <- length(TO_TEST.cov) + 1:length(OR.gt.9)
	TO_TEST <- 1:ncol(X.train)
	TO_INCL <- c()
	TO_OMIT <- c()
	N_PREDS <- ncol(X.train)
	MOD_SIZE <- min( MOD_SIZE, N_PREDS )
	AUC.test <- list()
	AUC.test$TS <- AUC.test$TR <- array( , c(N_PREDS,MOD_SIZE) )
	rownames(AUC.test$TS) <- rownames(AUC.test$TR) <- colnames(X.train)[1:N_PREDS]
	AUC <- array( , dim=c(MOD_SIZE,2) )
	colnames(AUC) <- c("TR","TS")
	# for ( r in 1:nrow(COEF.all) ) {
	print("## Building Models ##")
	start_time <- proc.time()
	for ( i in 1:MOD_SIZE ) {
		## Loop through all Predictors and Find Best on to ADD to model
		for ( r in TO_TEST ) {
			This_Pred <- colnames(X.train)[r]
			Which_Pred <- c( TO_INCL, This_Pred )
			## Specify Phenotype Values
			PHENO.TR.9 <- factor(Y.train[,1])
			PHENO.TS.9 <- factor(Y.test[,1])
			## Create Binary Matrix for Genotypes
			BIN.TR.1 <- data.matrix( X.train[, Which_GT ] ) ; BIN.TR.1[which(BIN.TR.1==2,arr.ind=T)] <- 1
			BIN.TS.1 <- data.matrix( X.test[, Which_GT ] ) ; BIN.TS.1[which(BIN.TS.1==2,arr.ind=T)] <- 1
			## Calculate Score for each Person
			SCOR.TR.1 <- BIN.TR.1 %*% COEF.all[Which_GT,"OR.all"]
			SCOR.TS.1 <- BIN.TS.1 %*% COEF.all[Which_GT,"OR.all"]


			## Build Temporary Data Frame with These Covariates
			TEMP.TR.9 <- data.matrix( X.train[, Which_Pred ] ) ; colnames(TEMP.TR.9) <- Which_Pred
			TEMP.TS.9 <- data.matrix( X.test[, Which_Pred ] ) ; colnames(TEMP.TS.9) <- paste("TEMP.TR.9",Which_Pred,sep="")
			## Fit Model Using lm()
			MOD.TR.9 <- glm( PHENO.TR.9 ~ TEMP.TR.9, family=binomial(link=logit) )
			COEF.mod <- coef(MOD.TR.9)
			PRED.TR.9 <- COEF.mod[1] + TEMP.TR.9 %*% matrix(COEF.mod[-1])
			PRED.TS.9 <- COEF.mod[1] + TEMP.TS.9 %*% matrix(COEF.mod[-1])
			## Calculate Probability from Link Function
			 # p <- exp( b0 + b1X1 + ... + bnXn ) / ( 1 + exp( b0 + b1X1 + ... + bnXn ) )
			PROB.TR.9 <- exp(PRED.TR.9) / ( 1 + exp(PRED.TR.9) )
			PROB.TS.9 <- exp(PRED.TS.9) / ( 1 + exp(PRED.TS.9) )


			## Build Temporary Data Frame with These Covariates
			if (GT_Flag==1) {
				TEMP.TR.1 <- data.matrix(data.frame( MG.TR.1[, c(Which_Cov) ], SCOR.TR.1 ))
				TEMP.TS.1 <- data.matrix(data.frame( MG.TS.1[, c(Which_Cov) ], SCOR.TS.1 ))
				}else{
				TEMP.TR.1 <- data.matrix( MG.TR.1[, c(Which_Cov) ] )
				TEMP.TS.1 <- data.matrix( MG.TS.1[, c(Which_Cov) ] )
			}
			## Fit Model Using lm()
			MOD.TR.1 <- glm( PHENO.TR.1 ~ TEMP.TR.1, family=binomial(link=logit) )
			COEF.mod <- coef(MOD.TR.1)
			PRED.TR.1 <- COEF.mod[1] + TEMP.TR.1 %*% matrix(COEF.mod[-1])
			PRED.TS.1 <- COEF.mod[1] + TEMP.TS.1 %*% matrix(COEF.mod[-1])
			## Calculate Probability from Link Function
			 # p <- exp( b0 + b1X1 + ... + bnXn ) / ( 1 + exp( b0 + b1X1 + ... + bnXn ) )
			PROB.TR.1 <- exp(PRED.TR.1) / ( 1 + exp(PRED.TR.1) )
			PROB.TS.1 <- exp(PRED.TS.1) / ( 1 + exp(PRED.TS.1) )
			## Loop through different thresholds for calling Case/Control
			THRESHOLDS.1 <- seq( THRSH_RANGE[1],THRSH_RANGE[2],.002 )
			THRESHOLDS.2 <- seq( THRSH_RANGE[1],THRSH_RANGE[2], length.out=20 )
			if (length(THRESHOLDS.1)>length(THRESHOLDS.2) ) {
				THRESHOLDS <- THRESHOLDS.1
			}else{ THRESHOLDS <- THRESHOLDS.2 }
			N_THRSH <- length(THRESHOLDS)
			TPR.TR.1 <- TPR.TS.1 <- FPR.TR.1 <- FPR.TS.1 <- numeric(N_THRSH)
			for ( t in 1:N_THRSH ) {
				thresh <- THRESHOLDS[t]
				## Predict Cases/Controls
				BIN.TR.1 <- PROB.TR.1 > thresh
				BIN.TS.1 <- PROB.TS.1 > thresh
				## Calculate False Positive and True Postive Rates
				TPR.TR.1[t] <- length(which( PHENO.TR.1==2 & BIN.TR.1==T )) / length(which( PHENO.TR.1==2 ))
				TPR.TS.1[t] <- length(which( PHENO.TS.1==2 & BIN.TS.1==T )) / length(which( PHENO.TS.1==2 ))
				FPR.TR.1[t] <- length(which( PHENO.TR.1==1 & BIN.TR.1==T )) / length(which( PHENO.TR.1==1 ))
				FPR.TS.1[t] <- length(which( PHENO.TS.1==1 & BIN.TS.1==T )) / length(which( PHENO.TS.1==1 ))
			}
			# plot( 0,0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="FPR",ylab="TPR", main="ROC" )
			# abline( 0,1, lty=1, col="black", lwd=2 )
			# abline( h=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
			# abline( v=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
			# points( FPR.TR.1, TPR.TR.1, type="o", pch=20, col="chartreuse2", lwd=2 )
			# points( FPR.TS.1, TPR.TS.1, type="o", pch=20, col="dodgerblue2", lwd=2 )
			## Calculate Area Under Curve
			AUC.test$TR[r,i] <- sum( -( FPR.TR.1[2:N_THRSH]-FPR.TR.1[2:N_THRSH-1] ) * ( TPR.TR.1[2:N_THRSH]+TPR.TR.1[2:N_THRSH-1] ) /2 )
			AUC.test$TS[r,i] <- sum( -( FPR.TS.1[2:N_THRSH]-FPR.TS.1[2:N_THRSH-1] ) * ( TPR.TS.1[2:N_THRSH]+TPR.TS.1[2:N_THRSH-1] ) /2 )
		} # Close "r" Loop
	} # Close "i" Loop
	CV_COMP[[iter]] <- AUC
} # Close "n" Loop


######################################################
## END OF DOC ########################################
######################################################








































# ## Make list of all SNPs
# ALL_VARS <- Reduce(union, list(AFR$ID, EUR$ID, ASN$ID, SAN$ID, AMR$ID, ALL$ID, TBL$ID) )
# TBL_VARS <- TBL$ID

## Make list compiling all other data
FULL <- list(ALL, EUR, AFR, ASN, SAN, AMR)
names(FULL) <- c("ALL","EUR","AFR","ASN","SAN","AMR")
for (d in 1:length(FULL)) {
	FULL[[d]][which(FULL[[d]]$POS==112241766),]$ID <- "rs671"
}

## Go thru Variants and Calculate "Pooled" Allele Frequencies (based on TBL collapsing)
n_vars <- nrow(TBL)
VAR_FREQ <- array(, dim=c(n_vars,16) )
rownames(VAR_FREQ) <- TBL$ID
colnames(VAR_FREQ) <- c("TBL_A","TBL_B","CEU_A","CEU_B","ALL_A","ALL_B","EUR_A","EUR_B","AFR_A","AFR_B","AMR_A","AMR_B","ASN_A","ASN_B","SAN_A","SAN_B")
for ( v in 1:n_vars ) {
	var <- as.character(TBL$ID[v])
	GT_A <- as.character(TBL$GT_A[v])
	GT_B <- as.character(TBL$GT_B[v])
	PH_A <- TBL$PHENO_A[v]
	PH_B <- TBL$PHENO_B[v]
	VAR_FREQ[v,"TBL_A"] <- TBL$TBL_FRQ_A[v]
	VAR_FREQ[v,"TBL_B"] <- TBL$TBL_FRQ_B[v]
	VAR_FREQ[v,"CEU_A"] <- TBL$CEU_FRQ_A[v]
	VAR_FREQ[v,"CEU_B"] <- TBL$CEU_FRQ_B[v]
	for ( d in 1:length(FULL) ) {
		DAT <- FULL[[d]]
		DAT_NM <- names(FULL)[d]
		if ( var %in% DAT$ID ) {
			vd <- which(DAT$ID==var)
			# Make Genotype Options
			HOM_REF <- paste(rep(DAT$REF[vd],2),collapse="")
			HET_1 <- paste(DAT$REF[vd],DAT$ALT[vd],sep="",collapse="")
			HET_2 <- paste(DAT$ALT[vd],DAT$REF[vd],sep="",collapse="")
			HOM_ALT <- paste(rep(DAT$ALT[vd],2),collapse="")
			# Check for presence of strand issue
			if ( grepl(HOM_REF,GT_A)==F & grepl(HOM_REF,GT_B)==F ) {
				HOM_REF <- STRAND(HOM_REF)
				HET_2 <- STRAND(HET_1)
				HET_1 <- STRAND(HET_2)
				HOM_ALT <- STRAND(HOM_ALT)
			}
			# Determine how to pool genotypes
			DAT_FQ_A <- 0
			DAT_FQ_B <- 0
			if ( grepl(HOM_REF, GT_A) ) { DAT_FQ_A <- DAT_FQ_A + DAT$HOM.REF[vd] }
			if ( grepl(HOM_REF, GT_B) ) { DAT_FQ_B <- DAT_FQ_B + DAT$HOM.REF[vd] }
			if ( grepl(HOM_ALT, GT_A) ) { DAT_FQ_A <- DAT_FQ_A + DAT$HOM.VAR[vd] }
			if ( grepl(HOM_ALT, GT_B) ) { DAT_FQ_B <- DAT_FQ_B + DAT$HOM.VAR[vd] }
			if ( grepl(HET_1, GT_A) | grepl(HET_2, GT_A) ) { DAT_FQ_A <- DAT_FQ_A + DAT$HET[vd] }
			if ( grepl(HET_1, GT_B) | grepl(HET_2, GT_B) ) { DAT_FQ_B <- DAT_FQ_B + DAT$HET[vd] }
			VAR_FREQ[v,paste(DAT_NM,"A",sep="_")] <- DAT_FQ_A
			VAR_FREQ[v,paste(DAT_NM,"B",sep="_")] <- DAT_FQ_B
		}
	}
}

for ( v in 1:n_vars ) {
	if ( all( is.na(VAR_FREQ[v,c(1:4,7:8)])==F ) ) {
		var <- as.character(TBL$ID[v])
		GT_A <- as.character(TBL$GT_A[v])
		GT_B <- as.character(TBL$GT_B[v])
		PH_A <- gsub("_"," ",TBL$PHENO_A[v])
		PH_B <- gsub("_"," ",TBL$PHENO_B[v])
		PH <- gsub("_"," ",TBL$PHENO[v])
		print(var)
		PNG_WRITE <- paste(DATE,"_",v,"_",var,".png",sep="")
		png(gsub("FILENAME",PNG_WRITE,PathToData), height=1200,width=1200,pointsize=27)
		PLOT_DAT <- round( rbind( VAR_FREQ[v,seq(1,ncol(VAR_FREQ),2)], VAR_FREQ[v,seq(2,ncol(VAR_FREQ),2)] ),0)
		barplot(PLOT_DAT[1,]/colSums(PLOT_DAT), col="chartreuse2", ylim=c(-1.2,1.2), las=2, width=.8,space=.25, main=paste("Allele Freq predicting", PH, "-",var), ylab=paste(PH_B,"<--",PH,"-->",PH_A) )
		barplot(-PLOT_DAT[2,]/colSums(PLOT_DAT), col="dodgerblue2", add=T, xaxt="n", yaxt="n", width=.8,space=.25 )
		arrows(.5,1.005,1.5,1.005, length=.1, code=3, angle=90, lwd=3)
		arrows(.5,1.105,3.5,1.105, length=.1, code=3, angle=90, lwd=3)
		# Run Test
		CHI_HAPMAP <- fisher.test(PLOT_DAT[,c(1,2)])$p.value
		if ( CHI_HAPMAP <= .05 & CHI_HAPMAP > .01 ) { text(1,1.025, labels="*") }
		if ( CHI_HAPMAP <= .01 & CHI_HAPMAP > .001 ) { text(1,1.025, labels="**") }
		if ( CHI_HAPMAP <= .001 ) { text(1,1.025, labels="***") }
		CHI_1KG_CEU <- fisher.test(PLOT_DAT[,c(1,4)])$p.value
		if ( CHI_1KG_CEU <= .05 & CHI_1KG_CEU > .01 ) { text(2,1.125, labels="*") }
		if ( CHI_1KG_CEU <= .01 & CHI_1KG_CEU > .001 ) { text(2,1.125, labels="**") }
		if ( CHI_1KG_CEU <= .001 ) { text(2,1.125, labels="***") }
		text(8,1,labels=v)
		text(5,1.1, labels=paste(GT_A,"-",PH_A))
		text(5,-1.1, labels=paste(GT_B,"-",PH_B))
		dev.off()
	}
}















######################################################
## END OF DOC ########################################
######################################################
