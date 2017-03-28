# 200K MCMC iterations W no ascert
library(mvtnorm)
library(parallel)
library(reshape)

library(MCMCglmm)
server <- strsplit(getwd(),"/")[[1]][2]
lib_folder <- paste0("/",server, "/mbatchou/GMMAT") 
gmmat_load <- require(GMMAT, lib.loc = lib_folder)
if (!gmmat_load) {
	lib_folder <- paste0("/",server, "/mbatchou/GMMAT_compute78")
	gmmat_load <- require(GMMAT, lib.loc = lib_folder)
}
if (!gmmat_load) library(GMMAT)

# Parameters for simulation
N_fam <- 45
n_person <- 22
mean_G1 <- 2 * 0.2
n_total <- N_fam * n_person

# Make pedigree
ped1 <- data.frame(fam = NA, ind = 1:n_total, fa = NA, mo = NA)
fa_ind <- c(NA,NA,NA,1,1,NA,1,NA,1,NA,4,4,4,6,6,6,8,8,10,10,10,10)
mo_ind <- c(NA,NA,NA,2,2,NA,2,NA,2,NA,3,3,3,5,5,5,7,7,9,9,9,9)
ped1[,c('mo','fa')] <- cbind(rep(mo_ind, N_fam)+ n_person*rep(0:(N_fam-1), each = n_person),
rep(fa_ind, N_fam)+ n_person*rep(0:(N_fam-1), each = n_person))
Prec <- inverseA(ped1[,c('ind','mo','fa')])$Ainv
ped1[is.na(ped1)] <- 0
ped1$fam <- rep(1:N_fam, each = n_person)
write.table(ped1,'pedigree',row.names=FALSE, col.names=FALSE,sep='\t')

# Generate data
phi <- as.matrix(read.table("add_phi_mat_n22")); colnames(phi) <- 1:n_person
chol_phi <- chol(phi); C_inv <- solve(chol_phi)
C_inv_t_1 <- colSums(C_inv)
H_mat <- diag(n_person) - matrix(rep(tcrossprod(C_inv_t_1, C_inv),n_person), n_person, n_person, byrow = T)/sum(C_inv_t_1^2)
Phimat <- diag(N_fam) %x% phi
kin <- melt(diag(N_fam)%x%phi)[,c(2,1,3)]; kin = kin[kin[,1]<=kin[,2],]
# Get fam for row and col indID
kin$famrow <- ((kin[,1]-1)%/%n_person) + 1; kin$famcol <- ((kin[,2]-1)%/%n_person) + 1
kin <- kin[kin$famrow==kin$famcol,]
kin2 <- kin; kin2$v[kin2$v == 1] <- 0; kin2$v <- kin2$v/2
phi_inv <- solve(phi) #22x22 matrix
sex <- rep(c(1,2,2,1,2,1,2,1,2,1,2,1,1,2,1,2,1,2,2,1,2,1),N_fam) #constant

peddrop <- function (parentmat, MAF) { # parentmat is matrix of parental info
	nfound <- sum((parentmat[, 1] == 0) & (parentmat[, 2] == 0))
	npeople <- dim(parentmat)[1]
	haplos <- matrix(0, npeople, 2)
	haplos[(parentmat[, 1] == 0) & (parentmat[, 2] == 0),] <- matrix(1:(2 * nfound), nfound, 2, byrow=T)

	for (i in 1:npeople) {
		#/* When both parents are in the pedigree drop alleles down the pedigree */
		if ((parentmat[i, 1] != 0) & (parentmat[i, 2] != 0)) {
			#/* Randomly choose chrom. 0 or 1 to transmit to kid at each locus
		trchrom = sample(1:2, 2, replace = T)
		#// For GT haplotype
		haplos[i, 1] = haplos[parentmat[i, 1], trchrom[1]];
		haplos[i, 2] = haplos[parentmat[i, 2], trchrom[2]];
		} 
	}

	state <- as.numeric(runif(2 * nfound) < MAF) #1 if minor, 0 if major allele
	genos <- matrix(0, npeople, 2)
	for (i in 1:npeople) genos[i, ] = state[haplos[i,]]

	return(genos)
}

Tu <- .2; Tot <- 6
add_var <- Tu * Tot
beta <- c(.05, 1, .5)
old_v_xb <- sum(c(363.6669,.25,1)*beta^2)
S_XB <- sqrt(Tot * (1 - Tu) / old_v_xb) # Ratio of target variance to var(X)
beta <- beta * S_XB
LG <- log(4) # increase of 10% in the odds of being a case
int <- qlogis(0.025) - beta[1] * 31.86 -  beta[2] * 1.5 - LG * mean_G1  # Prevalence - b * E(X) - lambda * mean of G1

runme <- function(dum){
	Y <- age <- Z <- G1 <- U_full <- NULL
# Generate family by family
	for (i_fam in 1:N_fam){
		ascert = 0
		while(ascert == 0){
			age_f <- c(73,75,46,43,40,46,40,43,47,51,18,21,15,15,12,9,13,17,24,21,18,14) + runif(n_person, -1.5, 1.5)
			Z_f <- rnorm(n_person)
			U <- c(rmvnorm(1, sigma = add_var * phi)) 
			G1_f <- rowSums(peddrop( ped1[1:n_person, c("fa","mo")], 0.2 )) # Assumes same pedigree structure across families
			
			Y_f <- rbinom(n_person, size = 1, prob = plogis(int + U + LG * G1_f + age_f * beta[1] + sex[ped1$fam == i_fam] * beta[2] + Z_f * beta[3]))
			
			if(sum(Y_f) >= 0) ascert <- 1 #No asc.
		}
		Y <- append(Y, Y_f)
		age <- append(age, age_f)
		Z <- append(Z, Z_f)
		U_full <- append(U_full, U)
		G1 <- append(G1, G1_f)
	}
	fam_data = data.frame(Y=Y, age = age, sex = sex, Z = Z, U = U_full, idU = 1:n_total)

### MCMC glmm
	prior1 <- list(R = list(V = 1, fix=1),
	G = list(G1 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=10^2)))
	mMCMC <- MCMCglmm(Y ~ 1 + age + sex + Z, random =~idU, family = "categorical", ginverse = list(idU = Prec), prior=prior1,data = fam_data, verbose = F, nitt=300e3, thin=500, burnin=15e3, pr = TRUE, pl = TRUE)
	U_hat_mode <- posterior.mode(mMCMC$Sol[,-(1:mMCMC$Fixed$nfl)])
	U_hat_mean <- colMeans(mMCMC$Sol[,-(1:mMCMC$Fixed$nfl)])
	ue <- as.matrix(mMCMC$Liab) - tcrossprod(as.matrix(mMCMC$Sol[,1:mMCMC$Fixed$nfl]), as.matrix(mMCMC$X))
	ue_mode <- posterior.mode(ue); ue_mean <- colMeans(ue)

## Assessment plots
	pdf(paste0("corGU",dum,".pdf"), width = 8.5, height = 6.5)
	par(mfrow=c(2,2))
	df_u <- crossprod(diag(N_fam)%x%C_inv,as.matrix(data.frame(U_hat_mean, U_hat_mode, ue_mean, ue_mode)))
	v <- sapply(1:ncol(df_u), function(i){
		boxplot(df_u[,i]~ G1)
		title(c("U_hat_mean", "U_hat_mode", "ue_mean", "ue_mode")[i])
	})
# ggplot(data = data.frame(U=fam_data$U, U_s, G1), aes(x=U_s, group = G1, fill = G1))+geom_histogram(position="identity",alpha=0.5,binwidth=0.1)+theme_bw()

	df_u <- crossprod(diag(N_fam)%x%C_inv, as.matrix(data.frame(U_hat_mean, U_hat_mode, ue_mean, ue_mode)) - fam_data$U)
	v <- sapply(1:ncol(df_u), function(i){
		boxplot(df_u[,i]~ G1)
		title(c("U_hat_mean", "U_hat_mode", "ue_mean", "ue_mode")[i])
	})
	
	df_u <- crossprod(diag(N_fam)%x%C_inv,as.matrix(data.frame(U_hat_mean, U_hat_mode, ue_mean, ue_mode)))
	v <- sapply(1:ncol(df_u), function(i){
		boxplot(df_u[,i]~ round(crossprod(diag(N_fam)%x%C_inv,G1),4))
		title(c("U_hat_mean", "U_hat_mode", "ue_mean", "ue_mode")[i])
	})
	dev.off()
	return(0)
}

mclapply(1:10, runme, mc.cores = 12)