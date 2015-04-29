## Compare Biggest Loser Genotype Frequencies w/ 1KG Data ##
## Filter Data by Relationships/Ancestry ##
## February 5, 2015 ##
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
DATE <- "20150428"
PathToPlot <- paste("/Users/kstandis/Data/TBL/Plots/",DATE,"_Plots",sep="")
PathToWrite <- "/Users/kstandis/Data/TBL/Data/Filtered_Tables"

## Path to Old Data (8/25)
PathToOldData <- "/Users/kstandis/Data/TBL/Data/FILENAME"

## Path to New Set of Data
PathToNewData <- "/Users/kstandis/Data/TBL/Data/20141001/FILENAME"

## Path to New Set of 1KG VCF files
PathTo1KG <- "/Users/kstandis/Data/TBL/Data/20141001/1KG/FILENAME"

## Path to List of SNPs w/ Strand Issues
PathToStrandSNPs <- "/Users/kstandis/Data/TBL/Data/20150323_Strand_Issue_SNPs.xlsx"

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
proc.time()-start

## Load rsID/Gene/Phenotype Info (10/01)
PGR_DIET <- read.xlsx(gsub("FILENAME","FIT genes.xlsx",PathToNewData), sheetName="Diet", colIndex=3:6, rowIndex=2:69, as.data.frame=T, header=T)
PGR_NUTR <- read.xlsx(gsub("FILENAME","FIT genes.xlsx",PathToNewData), sheetName="Nutrition", colIndex=3:6, rowIndex=2:11, as.data.frame=T, header=T)
PGR_EXER <- read.xlsx(gsub("FILENAME","FIT genes.xlsx",PathToNewData), sheetName="Exercise", colIndex=3:6, rowIndex=3:14, as.data.frame=T, header=T)
PGR_META <- read.xlsx(gsub("FILENAME","FIT genes.xlsx",PathToNewData), sheetName="Metabolic Health", colIndex=3:6, rowIndex=2:55, as.data.frame=T, header=T)
proc.time()-start

## Load Phenotype/Genotype/Clinical Data for Patients (10/01)
PHENO <- read.table(gsub("FILENAME","Phenotypic-Outcomes-Format-2.txt",PathToNewData), header=T, sep="\t" ) ; proc.time()-start
GENO <- read.table(gsub("FILENAME","Genotypes-Format-2.txt",PathToNewData), header=T, sep="\t", fill=T, comment.char="" ) ; proc.time()-start
CLIN.l <- read.xlsx(gsub("FILENAME","BL all sample list_130426_FINAL_pg1.xlsx",PathToNewData), sheetIndex=1, colIndex=1:139, rowIndex=1:289, as.data.frame=T, header=T)
CLIN_KEY <- data.frame( NM_1=colnames(CLIN.l), NM_2=c(CLIN.l[1,],recursive=T), row.names=NULL )
CLIN <- CLIN.l[2:nrow(CLIN.l),]
POOL <- read.table(gsub("FILENAME","SNP_Pheno_Pool_Summary.txt",PathToNewData), header=T, sep="\t", fill=T )
proc.time()-start

## Load 1KG files (10/01) 
 # Pedigree Info
PED <- read.table(gsub("FILENAME","1KG_Pedigree_20130606_g1k.ped",PathTo1KG), header=T, sep="\t")
INFO_PED <- read.xlsx(gsub("FILENAME","1KG_20130606_sample_info.xlsx",PathTo1KG), sheetName="Sample Info", colIndex=1:15, rowIndex=1:3501, as.data.frame=T, header=T)
proc.time()-start
 # VCF Files
VCF_ALL <- read.table(gsub("FILENAME","1KG_All_Vars.vcf",PathTo1KG), sep="\t", header=T, skip=254,comment.char="")
VCF_AFR <- read.table(gsub("FILENAME","1KG_AFR.recode.vcf",PathTo1KG), sep="\t", header=T, skip=254,comment.char="")
VCF_AMR <- read.table(gsub("FILENAME","1KG_AMR.recode.vcf",PathTo1KG), sep="\t", header=T, skip=254,comment.char="")
VCF_ASN <- read.table(gsub("FILENAME","1KG_ASN.recode.vcf",PathTo1KG), sep="\t", header=T, skip=254,comment.char="")
VCF_SAN <- read.table(gsub("FILENAME","1KG_SAN.recode.vcf",PathTo1KG), sep="\t", header=T, skip=254,comment.char="")
VCF_EUR <- read.table(gsub("FILENAME","1KG_EUR.recode.vcf",PathTo1KG), sep="\t", header=T, skip=254,comment.char="")
proc.time()-start
 # Ancestry Key
KEY <- read.table(gsub("FILENAME","Panel_Key.txt",PathTo1KG), header=T)

## Load List of Strand Issues
STRAND_LIST <- read.xlsx( PathToStrandSNPs, sheetIndex=1,rowIndex=1:19,colIndex=1:3,header=T )
colnames(STRAND_LIST) <- c("ID","Flip","Flip_2")
FLIP.SNPs <- as.character( STRAND_LIST$ID[ which(STRAND_LIST$Flip=="No") ] )
FLIP.SNPs <- gsub( " ","", FLIP.SNPs )

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
TROUBLE_SNPS[ which( TROUBLE_SNPS %in% FLIP.SNPs ) ]

## Plot correlation b/n Allele Frequency
png( paste(PathToPlot,"Scatter_Allele_Freq.png",sep="/"), height=1200,width=1200, pointsize=30)
plot( 0,0, type="n", xlab="1KG Allele Freq",ylab="TBL Allele Freq", main="MAF - TBL vs 1KG",xlim=c(0,1),ylim=c(0,1))
abline( h=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
abline( v=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
abline( 0,1, lty=2, col="black", lwd=2 )
COLS <- rep("dodgerblue2",nrow(POP)) ; COLS[TROUBLE] <- "firebrick2"
points(POP$AF_1KG,POP$AF_TBL, col=COLS, pch=20 )
legend( 0,1, col=c("dodgerblue2","firebrick2"), pch=20, legend=c("Confirmed","Potential Strand Issue") )
dev.off()

######################################################
## WRITE TABLES ######################################
######################################################

## Write Tables Needed for Analyses
write.table( POP_TBL_GC, paste(PathToWrite,"POP_TBL_GC.txt",sep="/"), sep="\t",row.names=F,col.names=T,quote=F )
write.table( POP_1KG_GC, paste(PathToWrite,"POP_1KG_GC.txt",sep="/"), sep="\t",row.names=F,col.names=T,quote=F )
write.table( VARS_BOTH, paste(PathToWrite,"VARS_BOTH.txt",sep="/"), sep="\t",row.names=F,col.names=F,quote=F )
write.table( POP, paste(PathToWrite,"POP.txt",sep="/"), sep="\t",row.names=T,col.names=T,quote=F )
write.table( TROUBLE_SNPS, paste(PathToWrite,"TROUBLE_SNPS.txt",sep="/"), sep="\t",row.names=F,col.names=F,quote=F )
write.table( POOL, paste(PathToWrite,"POOL.txt",sep="/"), sep="\t",row.names=F,col.names=T,quote=F )





































######################################################
## END OF DOC ########################################
######################################################
