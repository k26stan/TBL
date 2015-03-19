## Compare Biggest Loser Genotype Frequencies w/ 1KG Data ##
## Run single-locus association tests ##
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
DATE <- "20150205"
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

POP_TBL_GC <- 
POP_1KG_GC <- 
VARS_BOTH <- 
POP <- 






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
## END OF DOC ########################################
######################################################
