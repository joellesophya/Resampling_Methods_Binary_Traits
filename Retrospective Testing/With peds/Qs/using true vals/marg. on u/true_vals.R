# Would using marginal residual give better power??? -- we use true values for beta and sigma^2
# W ASCERT 6+
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
n_total <- N_fam * n_person
maf_G2 <- .2
mean_G1G2 <- (1 - .9^2) * (1-(1 - maf_G2)^2)

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

runme = function (dummy, Tu, Tot)
{
	kt<<- kt + 1; write(kt, "iteration_track.txt")
	add_var <- Tu * Tot
	betas <- c(.05, 1, .5)
	old_v_xb <- sum(c(363.6669,.25,1) * betas^2)
	S_XB <- sqrt(Tot * (1 - Tu) / old_v_xb) # Ratio of target variance to var(X)
	betas <- betas * S_XB
	LG <- log(1.6) # increase of 10% in the odds of being a case
	int <- qlogis(0.025) - betas[1] * 31.86 -  betas[2] * 1.5 - LG * 2 * .1  # Prevalence - b * E(X) - lambda * mean of epistatis indicator
	pvals_t1err = NA; pvals_power = NA
	
	while(any(is.na(pvals_t1err)) | any(is.na(pvals_power))){
		Y <- age <- Z <- G1 <- u <- NULL  
		# G2 <- NULL  
		# Generate family by family
		for (i_fam in 1:N_fam){
			ascert = 0
			while(ascert == 0){
				age_f <- c(73,75,46,43,40,46,40,43,47,51,18,21,15,15,12,9,13,17,24,21,18,14) + runif(n_person, -1.5, 1.5)
				Z_f <- rnorm(n_person)
				U <- c(rmvnorm(1, sigma = add_var * phi)) 
				G1_f <- rowSums(peddrop( ped1[1:n_person, c("fa","mo")], 0.1 )) # Assumes same pedigree structure across families
				# G2_f <- rowSums(peddrop( ped1[1:n_person, c("fa","mo")], maf_G2 ))
				
				Y_f <- rbinom(n_person, size = 1, prob = plogis(int + U + LG * G1_f +	age_f * betas[1] + sex[ped1$fam == i_fam] * betas[2] + Z_f * betas[3]))
				
				if(sum(Y_f) > 5) ascert <- 1 #At least 5 affected
			}
			Y <- append(Y, Y_f)
			age <- append(age, age_f)
			Z <- append(Z, Z_f)
			u <- append(u, U)
			G1 <- append(G1, G1_f)
			# G2 <- append(G2, G2_f)
		}
		fam_data = data.frame(Y=Y, age = age, sex = sex, Z = Z, idU = 1:n_total)
		
		Ygu <- LG * G1 + u
		target_var <- summary(lm(crossprod(diag(N_fam) %x% C_inv,Ygu)~crossprod(diag(N_fam) %x% C_inv,rep(1,n_total))-1))$sig^2

		### MCMC glmm
		ck <- ((16*sqrt(3))/(15*pi))^2
		Xb <- as.vector(model.matrix(~age + sex + Z, data = fam_data) %*% c(int,betas))
		pi_hat <- plogis(Xb / sqrt(1 + ck * target_var))

## Generating unassociated Gs
		n_G <-50 # 50 * n_rep replicates for type 1 error
		Gmat <-  cbind(G1, replicate(n_G, rowSums(peddrop( ped1[,c("fa","mo")], 0.3 )))) # G1 is first column + n * n_G matrix 
# Numerator
		res_fam <- matrix(Y - pi_hat, n_person, N_fam)
		T_num_G <- colSums(Gmat * c(crossprod(H_mat, res_fam)))^2
		
# Retrospective denominator
		get_denum_retro <- function(G, Matdenum){ 
			mafs <- colMeans(G) / 2
			2 * mafs * (1 - mafs) * Matdenum
		}
		YHPhiY <- sum(c(crossprod(H_mat, res_fam)) * c(crossprod(phi, res_fam)))
		T_denum_G <- get_denum_retro(Gmat, YHPhiY)
		pvalstot <- pchisq(T_num_G / T_denum_G, df = 1, lower.tail = F)
		pvals_t1err <- pvalstot[-1]
		pvals_power <- pvalstot[1]

		numend <- paste0("_",sample(1e6,1),"_")
		filend <- paste0(numend, '.txt')
# GMMAT
		mGMMAT <- glmmkin(Y ~ age + sex + Z,  data = fam_data, kins = Phimat, family = binomial(link="logit"))
		write.table(cbind(paste0("SNP", 1:(n_G+1)), t(Gmat)), paste0('genofile', filend), row.names=F, col.names=F, quote=F, sep='\t')
		capture.output(glmm.score(mGMMAT, infile = paste0('genofile', filend), outfile = paste0('gmmat_out', filend)), file = "/dev/null")
		if(file.exists(paste0('gmmat_out', filend))){
		pvalstot <- as.numeric(read.table(paste0('gmmat_out', filend), header=T)[,6])
		pvals_t1err <- cbind(pvals_t1err, pvalstot[-1]) 
		pvals_power <- append(pvals_power, pvalstot[1]) 
		} else{
		pvals_t1err <- cbind(pvals_t1err, NA) 
		pvals_power <- append(pvals_power, NA) 
		write(paste('Error in GMMAT', kt), "errorfile.txt", append = T)}
		system(paste0('rm -R *', numend,'*'))

# CERAMIC
		system(paste('mkdir', numend))	
		write.table(cbind(ped1[,1:2],fam_data$Y+1,fam_data[,2:4]), paste0('./',numend,'/','phenofile'), sep = '\t', quote=F, col.names=F, row.names=F)
		kin2 = kin; kin2$value[kin2$v == 1] = 0; kin2$value = kin2$v/2
		write.table(kin2[,c(4,1:3)],paste0('./',numend,'/','grmfile'), row.names=F, col.names=F,sep=' ')
		write(Gmat, paste0('./',numend,'/','genofile'), ncolumns = n_total, sep='\t')
		com_ceramic <- paste0('(cd ./',numend,' && exec ../BATMAN -p phenofile -g genofile -k grmfile -c)')
		ret <- system(com_ceramic, ignore.stdout = T, ignore.stderr = T)
		if ((!file.exists(paste0('./',numend,'/BATMANtest.pvalues'))) | (ret !=0)){
			pvals_t1err <- cbind(pvals_t1err, NA) 
			pvals_power <- append(pvals_power, NA) 
			write(paste('Error in CERAMIC', kt), "errorfile.txt", append = T)
		} else{
			pvalstot <- as.numeric(read.table(paste0('./',numend,'/BATMANtest.pvalues'), header=T)[,3])
			pvals_t1err <- cbind(pvals_t1err, pvalstot[-1]) 
			pvals_power <- append(pvals_power, pvalstot[1]) 
		}
		system(paste0('rm -R *', numend,'*'))
	}
	pvalstot <- list(Type1Err = pvals_t1err, Power = pvals_power)
    return(pvalstot)
}

M = rbind(c(0.2,5.5),c(0.4,6),c(.8,6.3),c(1,6.4))
n_M = nrow(M)
# Will be a list (over Tu/Te) of a list (over rep) of matrices (GMMAT and MCMC)
res = rep(list(list()), n_M) 
n_rep = 750; kt <- 0 # To keep track of iterations

for(i in 1:n_M) {
  res[[i]] = mclapply(1:n_rep, runme, mc.cores = 15, Tu= M[i,1], Tot = M[i,2])
  save.image() #Save the res object
}

allfile <- list.files()
file.remove(allfile[grep("BATMAN.",allfile)])
  
body <- tail(strsplit(getwd(),"/")[[1]],1); body <- paste("Folder: ", body) 
system(paste("echo \"",body,"\" | mailx -s \"Sim done\" joelle.mbatchou@gmail.com"), ignore.stdout = T, ignore.stderr = T)