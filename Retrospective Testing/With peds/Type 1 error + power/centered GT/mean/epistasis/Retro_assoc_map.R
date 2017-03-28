# Corrected denominator retro. test
# Include GMMAT analysis
# 250K MCMC iterations w no asc.
# GLS estimator for sigma_G in denum of retro. test
library(mvtnorm)
library(parallel)
library(reshape)
library(MCMCglmm)
library(data.table)
server <- strsplit(getwd(),"/")[[1]][2]
lib_folder <- paste0("/",server, "/mbatchou/GMMAT") 
gmmat_load <- require(GMMAT, lib.loc = lib_folder)
if (!gmmat_load) {
	lib_folder <- paste0("/",server, "/mbatchou/GMMAT_compute78")
	gmmat_load <- require(GMMAT, lib.loc = lib_folder)
}
if (!gmmat_load) library(GMMAT)
source("predict_func_mod.R")

# Parameters for simulation
N_fam <- 45
n_person <- 22
n_total <- N_fam * n_person
n_asc <- 4 # Degree of ascertainment
ncov <- 4
maf_G2 <- .2
mean_G1G2 <- (1 - .9^2) * (1-(1 - maf_G2)^2)


# Make pedigree
ped1 <- data.table(fam = NA, ind = 1:n_total, fa = NA, mo = NA)
fa_ind <- c(NA,NA,NA,1,1,NA,1,NA,1,NA,4,4,4,6,6,6,8,8,10,10,10,10)
mo_ind <- c(NA,NA,NA,2,2,NA,2,NA,2,NA,3,3,3,5,5,5,7,7,9,9,9,9)
ped1[,`:=`(mo = (rep(mo_ind, N_fam) + n_person*rep(0:(N_fam-1), each = n_person)),
fa = (rep(fa_ind, N_fam)+ n_person*rep(0:(N_fam-1), each = n_person)), fam = rep(1:N_fam, each = n_person))]
Prec <- inverseA(ped1[,.(ind,mo,fa)])$Ainv
ped1[is.na(ped1)] <- 0
write.table(ped1,'pedigree',row.names=FALSE, col.names=FALSE,sep='\t')
sex <- rep(c(1,2,2,1,2,1,2,1,2,1,2,1,1,2,1,2,1,2,2,1,2,1),N_fam) #constant

# Generate data
phi <- as.matrix(read.table("add_phi_mat_n22")); colnames(phi) <- 1:n_person
Phimat <- diag(N_fam)%x%phi
Phi_inv <- diag(N_fam)%x% chol2inv(chol(phi))
P <- Phi_inv - outer(rowSums(Phi_inv), colSums(Phi_inv)) / sum(Phi_inv)
setnames(kin <- data.table(melt(Phimat)), 1:3, c("person2","person1","kincoef"))
# setcolorder(kin, c("person1","person2","kin"))
kin <- kin[person1<=person2 & ((person1-1)%/%n_person) == ((person2-1)%/%n_person)]
kin[, fam := (person1-1)%/%n_person+ 1]
# For CERAMIC -- same for all realizations
kin[kincoef == 1, kincoef:= 0]; kin[, kincoef := kincoef/2]
write.table(kin[,.(fam, person1, person2, kincoef)],'grmfile', row.names=F, col.names=F,sep=' ')

peddrop <- function (parentmat, MAF) { # parentmat is matrix of parental info
	nfound <- parentmat[fa == 0 & mo == 0, .N]
	npeople <- nrow(parentmat)
	haplos <- matrix(0, npeople, 2)
	haplos[(parentmat[, fa] == 0) & (parentmat[, mo] == 0),] <- matrix(1:(2 * nfound), nfound, 2, byrow=T)

	for (i in 1:npeople) {
		#/* When both parents are in the pedigree drop alleles down the pedigree */
		if (sum(parentmat[i]) != 0) {
			#/* Randomly choose chrom. 0 or 1 to transmit to kid at each locus
			trchrom = sample(1:2, 2, replace = T)
			#// For GT haplotype
			haplos[i, 1] <- haplos[parentmat[i, fa], trchrom[1]]
			haplos[i, 2] <- haplos[parentmat[i, mo], trchrom[2]]
		} 
	}

	state <- as.numeric(runif(2 * nfound) < MAF) #1 if minor, 0 if major allele
	genos <- state[haplos]
	genos <- genos[1:npeople] + genos[-(1:npeople)] 
	return(genos)
}

runme = function (dummy, Tu, Tot)
{
	write(dummy, "iteration_track.txt")
	add_var <- Tu * Tot
	betas <- c(.05, 1, .5)
	old_v_xb <- sum(c(363.6669,.25,1) * betas^2)
	S_XB <- sqrt(Tot * (1 - Tu) / old_v_xb) # Ratio of target variance to var(X)
	betas <- betas * S_XB
	LG <- log(4.8) # increase of 10% in the odds of being a case
	
	int <- qlogis(0.025) - sum(betas * c(31.86, 1.5, 0)) - LG * mean_G1G2 # Prevalence - b * E(X) - lambda * mean of epistatis indicator
	U <- c(t(rmvnorm(N_fam, sigma = add_var * phi)))
	age <- rep(c(73,75,46,43,40,46,40,43,47,51,18,21,15,15,12,9,13,17,24,21,18,14), N_fam) + runif(n_total, -1.5, 1.5)
	Z <- rnorm(n_total)
	G1 <- peddrop(ped1[, .(fa,mo)], 0.1) 
	G2 <- peddrop(ped1[, .(fa,mo)], maf_G2) 
	###--###
	covars <- int + U  +	cbind(age, sex, Z) %*% betas + LG * ((G1 > 0) & (G2 > 0))
	pvalstot <- NULL
	
	while(is.null(pvalstot)){
		Y <- c(sapply(1:N_fam, function(i_fam){
			Y_f <- -5
			while(sum(Y_f) < n_asc) Y_f <- rbinom(n_person, size = 1, prob = plogis(covars[1:n_person + (i_fam - 1) * n_person]))
			return(Y_f)
		}))
		
		fam_data <- data.table(ped1, sex, Y, age, Z, idU = 1:n_total)
		numend <- paste0(dummy,"_",sample(1e6,1))
		filend <- paste0(numend, '.txt')
	
### MCMC glmm
		prior1 <- list(R = list(V = 1, fix=1),
		G = list(G1 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=10^2)))
		mMCMC <- MCMCglmm(Y ~ 1 + age + sex + Z, random =~idU, family = "categorical",
		ginverse = list(idU = Prec), prior=prior1,
		data = fam_data, verbose = F, nitt=250e3, thin=500, burnin=10e3, pr = TRUE)
		ck <- ((16*sqrt(3))/(15*pi))^2
		effS <- effectiveSize(mMCMC$VCV)[1]
		est_betas <- as.matrix(mMCMC$Sol[,1:ncov]) /sqrt(1 + ck)#niter*ncov
		est_u <- as.matrix(mMCMC$Sol[,-(1:ncov)]) /sqrt(1 + ck)
		est_s2 <- drop(mMCMC$VC[,1])/(1 + ck)
		
		## get post. mean of mu_y
		M1 <- Y - colMeans(plogis(as.matrix(tcrossprod(est_betas, mMCMC$X) + est_u)))
		M1 <- M1 - mean(M1)
		## get post mean of (b,u) to get mu_y
		M2 <- Y - plogis(drop(mMCMC$X %*% colMeans(est_betas) + colMeans(est_u)))
		M2 <- M2 - mean(M2)
		
## Generating unassociated Gs
		n_G <- 100 # m * n_rep replicates for type 1 error
		Gmat <-  cbind(G1, replicate(n_G, peddrop(ped1[, .(fa,mo)], 0.3))) # G1 is first column + n * n_G matrix 
		mafs <- colMeans(Gmat) / 2

# Numerator
		N1<- drop(crossprod(Gmat, M1)^2)
		N2<- drop(crossprod(Gmat, M2)^2)
		
# Retrospective denominator
		GLS_sig <- diag(crossprod(Gmat, P) %*% Gmat)/(n_total - 1)
		D1 <- GLS_sig * (M1 %*% crossprod(Phimat, M1))
		D2 <- GLS_sig * (M2 %*% crossprod(Phimat, M2))

		pvalstot <- pchisq(N1 / D1, df = 1, lower.tail = F)
		pvalstot <- cbind(pvalstot, pchisq(N2 / D2, df = 1, lower.tail = F))

# GMMAT
		mGMMAT <- glmmkin(Y ~ age + sex + Z,  data = fam_data, kins = Phimat, family = binomial(link="logit"))
		write.table(cbind(paste0("SNP", 1:(n_G+1)), t(Gmat)), paste0('genofile', filend), row.names=F, col.names=F, quote=F, sep='\t')
		capture.output(glmm.score(mGMMAT, infile = paste0('genofile', filend), outfile = paste0('gmmat_out', filend)), file = "/dev/null")
		if(file.exists(paste0('gmmat_out', filend))){
		pvalsG <- as.numeric(read.table(paste0('gmmat_out', filend), header=T)[,6])
		pvalstot <- cbind(pvalstot, pvalsG) 
		} else{
		pvalstot <- NULL
		write(paste('Error in GMMAT', dummy), "errorfile.txt", append = T)
		}

# CERAMIC
		system(paste('mkdir', numend))	
		write.table( fam_data[,.(fam, ind, Y+1, age, sex, Z)], paste0('./',numend,'/','phenofile'), sep = '\t', quote=F, col.names=F, row.names=F)
		write(Gmat, paste0('./',numend,'/','genofile'), ncolumns = n_total, sep='\t')
		com_ceramic <- paste0('(cd ./',numend,' && exec ../BATMAN -p phenofile -g genofile -k ../grmfile -c)')
		ret <- system(com_ceramic, ignore.stdout = T, ignore.stderr = T)
		if ((!file.exists(paste0('./',numend,'/BATMANtest.pvalues'))) | (ret !=0)){
			pvalstot <- NULL
			write(paste('Error in CERAMIC', dummy), "errorfile.txt", append = T)
		} else{
			pvalsC <- as.numeric(read.table(paste0('./',numend,'/BATMANtest.pvalues'), header=T)[,3])
			pvalstot <- cbind(pvalstot, pvalsC) 
		}
		system(paste0('rm -R *', numend, '*'))
# Erase files
	}
	pvalstot <- list(pvals = pvalstot, effectiveSize = effS)
	return(pvalstot)
}

M = data.frame(Tu = seq(0,1,.2), Tot = c(5.2,5.5,6,6.1,6.3,6.4))
n_M = nrow(M)
# Will be a list (over Tu/Te) of a list (over rep) of matrices (GMMAT and MCMC)
res = rep(list(list()), n_M) 
n_rep = 250

for(i in 1:n_M) {
	res[[i]] = mclapply(1:n_rep, runme, mc.cores = 16,  Tu= M$Tu[i], Tot = M$Tot[i])
	save.image() #Save the res object
}

allfile <- list.files()
file.remove(allfile[grep("BATMAN.",allfile)])

body <- tail(strsplit(getwd(),"/")[[1]],1); body <- paste("Folder: ", body) 
system(paste("echo \"",body,"\" | mailx -s \"Sim done\" joelle.mbatchou@gmail.com"), ignore.stdout = T, ignore.stderr = T)