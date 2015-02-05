## Compare Biggest Loser Genotype Frequencies w/ 1KG Data ##
## Do Multi-Locus/Polygenic Association Tests ##
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

# for ( n in 1:N_ITER) {
# 	AUC.plot <- CV_COMP[[n]]
# 	XLIM <- c(1,nrow(AUC.plot))
# 	YLIM <- c(0,1)
# 	png( paste(PathToOut,"/Poly_Train_",n,"_9_AUC.png",sep=""), height=1200,width=1600, pointsize=24)
# 	plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, main="AUC vs Model", ylab="AUC", xlab="",xaxt="n")
# 	axis( 1, at=1:nrow(AUC.plot), label=rownames(AUC.plot), las=2)
# 	abline( h=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
# 	abline( h=c(0,.5,1), lty=1, col="grey20", lwd=1 )
# 	abline( v=seq(0,XLIM[2],5), lty=2, col="grey50", lwd=1 )
# 	points( 1:nrow(AUC.plot), AUC.plot[,1], pch=20,type="o", col="chartreuse2", lwd=3)
# 	points( 1:nrow(AUC.plot), AUC.plot[,2], pch=20,type="o", col="deepskyblue2", lwd=3)
# 	dev.off()
# }

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
