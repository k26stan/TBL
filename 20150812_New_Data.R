## Check out New Data for TBL Data ##
## August 12, 2015 ##
## Kristopher Standish ##

#####################################################
## LOAD DATA ########################################
#####################################################

## Load Original Data & Re-Write in TXT files
# TAB.1 <- read.xlsx( "Data/TBL/Data/20150812/BL 130.xlsx", sheetIndex=1 )
# write.table( TAB.1, file="Data/TBL/Data/20150812/BL_130.1.txt",sep="\t",row.names=F,col.names=T,quote=F )
# TAB.2 <- read.xlsx( "Data/TBL/Data/20150812/BL 130.xlsx", sheetIndex=2 )
# write.table( TAB.2, file="Data/TBL/Data/20150812/BL_130.2.txt",sep="\t",row.names=F,col.names=T,quote=F )
# TAB.3 <- read.xlsx( "Data/TBL/Data/20150812/BL 130.xlsx", sheetIndex=3 )
# write.table( TAB.3, file="Data/TBL/Data/20150812/BL_130.3.txt",sep="\t",row.names=F,col.names=T,quote=F )

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to Files
PathToData <- "/Users/kstandis/Data/TBL/Data/20150812/"
PathToPlot <- paste("/Users/kstandis/Data/TBL/Plots/",DATE,sep="")
dir.create(PathToPlot)

## Load Data Tables
TAB.1 <- read.table( paste( PathToData,"BL_130.1.txt",sep=""), sep="\t",header=T, stringsAsFactors=F )
TAB.2 <- read.table( paste( PathToData,"BL_130.2.txt",sep=""), sep="\t",header=T, stringsAsFactors=F )
TAB.3 <- read.table( paste( PathToData,"BL_130.3.txt",sep=""), sep="\t",header=T, stringsAsFactors=F )

#####################################################
## FILTER/MERGE DATA TABLES #########################
#####################################################

## Filter out Empty Columns
TAB.1 <- TAB.1[ , -ncol(TAB.1) ]
TAB.2 <- TAB.2[ , -c(13,14) ]

## Merge Tables 2 & 3
TAB.23 <- merge( TAB.2, TAB.3, by.x="AccessionNumber", by.y="Accession.Number")

#####################################################
## ORGANIZE DATA ####################################
#####################################################

## Standardize Height and Weight Measurements #######

## Height (to Inches)
HT.key <- c( 1, 2.54, .0254 ) ; names(HT.key) <- c("inch","centimeter","meter")
 # Table 1
HT.1.a <- TAB.1[,c("AccessionNumber","Height","Height.Units")]
HT.1.b <- HT.1.a[,"Height"] / HT.key[ HT.1.a[,"Height.Units"] ] ; names(HT.1.b) <- HT.1.a[,1]
plot( HT.1.a[,"Height"],HT.1.b, col=as.numeric(as.factor(HT.1.a[,"Height.Units"])) )
 # Table 2/3
HT.23.a <- TAB.23[,c("AccessionNumber","Height.x","Height.Units.x")]
HT.23.b <- HT.23.a[,"Height.x"] / HT.key[ HT.23.a[,"Height.Units.x"] ] ; names(HT.23.b) <- HT.23.a[,1]
 # Remove Missing Height Values
HT.1.rm <- which(is.na(HT.1.b))
HT.1 <- HT.1.b[ -HT.1.rm ]
HT.23.rm <- which(is.na(HT.23.b))
HT.23 <- HT.23.b[ -HT.23.rm ]

## Weight ()
WT.key <- c( 1, 0.453592 ) ; names(WT.key) <- c("pound","kilogram")
 # Table 1
WT.1.a <- TAB.1[,c("AccessionNumber","Weight","Weight.Units")]
WT.1.b <- WT.1.a[,"Weight"] / WT.key[ WT.1.a[,"Weight.Units"] ] ; names(WT.1.b) <- WT.1.a[,1]
plot( WT.1.a[,"Weight"],WT.1.b, col=as.numeric(as.factor(WT.1.a[,"Weight.Units"])) )
 # Table 2/3
WT.23.a <- TAB.23[,c("AccessionNumber","Weight","Weight.Units.x")]
WT.23.b <- WT.23.a[,"Weight"] / WT.key[ WT.23.a[,"Weight.Units.x"] ] ; names(WT.23.b) <- WT.23.a[,1]
 # Remove Missing Height Values
WT.1.rm <- which(is.na(WT.1.b))
WT.1 <- WT.1.b[ -WT.1.rm ]
WT.23.rm <- which(is.na(WT.23.b))
WT.23 <- WT.23.b[ -WT.23.rm ]

## Combine Standardized HT/WT Data
SIZE.1 <- merge( data.frame(HT.1), data.frame(WT.1), by="row.names" )
BMI.1 <- WT.key["kilogram"]*SIZE.1$WT.1 / (HT.key["meter"]*SIZE.1$HT.1)^2
SIZE.1 <- data.frame( SIZE.1, BMI.1 ) ; colnames(SIZE.1) <- c("AccessionNumber","HT","WT","BMI")
SIZE.23 <- merge( data.frame(HT.23), data.frame(WT.23), by="row.names" )
BMI.23 <- WT.key["kilogram"]*SIZE.23$WT.23 / (HT.key["meter"]*SIZE.23$HT.23)^2
SIZE.23 <- data.frame( SIZE.23, BMI.23 ) ; colnames(SIZE.23) <- c("AccessionNumber","HT","WT","BMI")
# TEMP.1 <- merge( SIZE.1, TAB.1[,c("AccessionNumber","Calculated.BMI")] )
# TEMP.23 <- merge( SIZE.23, TAB.23[,c("AccessionNumber","Calculated.BMI")] )
# plot( TEMP.1$BMI, TEMP.1$Calculated.BMI, xlim=c(0,100),ylim=c(0,100) )
# plot( TEMP.23$BMI, TEMP.23$Calculated.BMI )

## Filter out Insane Numbers
RM.1 <- sort(Reduce( union, apply( SIZE.1[,2:4], 2, function(x) which( abs(x) > mean(x)+7*sd(x) ) ) ))
if ( length(RM.1)>0 ) { SIZE.1 <- SIZE.1[-RM.1, ] }
RM.23 <- sort(Reduce( union, apply( SIZE.23[,2:4], 2, function(x) which( abs(x) > mean(x)+7*sd(x) ) ) ))
if ( length(RM.23)>0 ) { SIZE.23 <- SIZE.23[-RM.23, ] }

## Merge w/ Genotype Results ########################
TAB.1.colnames <- 1:6
DAT.1 <- merge( SIZE.1, TAB.1[,TAB.1.colnames], by="AccessionNumber" )
TAB.23.colnames <- 1:6
DAT.23 <- merge( SIZE.23, TAB.23[,TAB.23.colnames], by="AccessionNumber" )

#####################################################
## CHECK OUT GENOTYPES ##############################
#####################################################

MOD <- list()
## ANOVA ############################################

MOD$ANOVA <- list()
## TABLE 1 ##
 # BMI vs Genotypes
MOD$ANOVA$BMI <- list()
MOD$ANOVA$BMI$F1 <- lm( BMI ~ FTO.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$ANOVA$BMI$M1 <- lm( BMI ~ MC4R.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$ANOVA$BMI$B1 <- lm( BMI ~ MC4R.genotype+FTO.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$ANOVA$BMI$I1 <- lm( BMI ~ MC4R.genotype*FTO.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
boxplot( BMI ~ MC4R.genotype, data=DAT.1, subset=DAT.1$MC4R.genotype!="--" )
boxplot( BMI ~ FTO.genotype, data=DAT.1, subset=DAT.1$FTO.genotype!="--" )
boxplot( BMI ~ MC4R.genotype+FTO.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--"), col=1:3, las=2 )
 # WT vs Genotypes
MOD$ANOVA$WT <- list()
MOD$ANOVA$WT$F1 <- lm( WT ~ FTO.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$ANOVA$WT$M1 <- lm( WT ~ MC4R.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$ANOVA$WT$B1 <- lm( WT ~ MC4R.genotype+FTO.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$ANOVA$WT$I1 <- lm( WT ~ MC4R.genotype*FTO.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
boxplot( WT ~ MC4R.genotype, data=DAT.1, subset=DAT.1$MC4R.genotype!="--" )
boxplot( WT ~ FTO.genotype, data=DAT.1, subset=DAT.1$FTO.genotype!="--" )
boxplot( WT ~ +MC4R.genotype+FTO.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--"), col=1:3, las=2 )
 # HT vs Genotypes
MOD$ANOVA$HT <- list()
MOD$ANOVA$HT$F1 <- lm( HT ~ FTO.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$ANOVA$HT$M1 <- lm( HT ~ MC4R.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$ANOVA$HT$B1 <- lm( HT ~ MC4R.genotype+FTO.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$ANOVA$HT$I1 <- lm( HT ~ MC4R.genotype*FTO.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
boxplot( HT ~ MC4R.genotype, data=DAT.1, subset=DAT.1$MC4R.genotype!="--" )
boxplot( HT ~ FTO.genotype, data=DAT.1, subset=DAT.1$FTO.genotype!="--" )
boxplot( HT ~ MC4R.genotype+FTO.genotype, data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--"), col=1:3, las=2 )

lapply( MOD$ANOVA, function(x) lapply( x, function(y) anova(y) ) )

## LINEAR MODEL #####################################

MOD$LIN <- list()
## BMI vs Genotypes
MOD$LIN$BMI$M1 <- lm( BMI ~ as.numeric(as.factor(MC4R.genotype)), data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$LIN$BMI$F1 <- lm( BMI ~ as.numeric(as.factor(FTO.genotype)), data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$LIN$BMI$B1 <- lm( BMI ~ as.numeric(as.factor(FTO.genotype))+as.numeric(as.factor(MC4R.genotype)), data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$LIN$BMI$I1 <- lm( BMI ~ as.numeric(as.factor(FTO.genotype))*as.numeric(as.factor(MC4R.genotype)), data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
plot( BMI ~ as.numeric(as.factor(FTO.genotype)), data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
abline( MOD$LIN$BMI$F1 )
## WT vs Genotypes
MOD$LIN$WT$M1 <- lm( WT ~ as.numeric(as.factor(MC4R.genotype)), data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$LIN$WT$F1 <- lm( WT ~ as.numeric(as.factor(FTO.genotype)), data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$LIN$WT$B1 <- lm( WT ~ as.numeric(as.factor(FTO.genotype))+as.numeric(as.factor(MC4R.genotype)), data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$LIN$WT$I1 <- lm( WT ~ as.numeric(as.factor(FTO.genotype))*as.numeric(as.factor(MC4R.genotype)), data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
## HT vs Genotypes
MOD$LIN$HT$M1 <- lm( HT ~ as.numeric(as.factor(MC4R.genotype)), data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$LIN$HT$F1 <- lm( HT ~ as.numeric(as.factor(FTO.genotype)), data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$LIN$HT$B1 <- lm( HT ~ as.numeric(as.factor(FTO.genotype))+as.numeric(as.factor(MC4R.genotype)), data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )
MOD$LIN$HT$I1 <- lm( HT ~ as.numeric(as.factor(FTO.genotype))*as.numeric(as.factor(MC4R.genotype)), data=DAT.1, subset=which(DAT.1$FTO.genotype!="--" & DAT.1$MC4R.genotype!="--") )

lapply( MOD$LIN, function(x) lapply( x, function(y) anova(y) ) )





