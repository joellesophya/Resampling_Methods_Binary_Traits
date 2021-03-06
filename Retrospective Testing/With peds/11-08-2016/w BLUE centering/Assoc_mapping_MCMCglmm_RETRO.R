# 200K MCMC iterations w no ascertainment
library(mvtnorm)
library(parallel)
library(reshape)

library(MCMCglmm)
server <- strsplit(getwd(),"/")[[1]][2]
gmmat_load <- tryCatch(suppressWarnings(library(GMMAT)), error = function(dummy){
lib_folder <- paste0("/",server, "/mbatchou/GMMAT_compute78")
library(GMMAT, lib.loc = lib_folder)
return(0)}
)

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

# Generate data
phi <- as.matrix(read.table("add_phi_mat_n22")); colnames(phi) <- 1:n_person
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
	beta <- c(.05, 1, .5)
	old_v_xb <- sum(c(363.6669,.25,1)*beta^2)
	S_XB <- sqrt(Tot * (1 - Tu) / old_v_xb) # Ratio of target variance to var(X)
	beta <- beta * S_XB
	LG <- log(4.8) # increase in penetrance
	int <- qlogis(0.025) - beta[1] * 31.86 -  beta[2] * 1.5 - LG * mean_G1G2  # Prevalence - b * E(X) - lambda * mean of epistatis indicator
	
	Y <- NULL; age <- NULL; Z <- NULL; G1 <- NULL  
	# G2 <- NULL  
	# Generate family by family
	for (i_fam in 1:N_fam){
		ascert = 0
		while(ascert == 0){
			age_f <- c(73,75,46,43,40,46,40,43,47,51,18,21,15,15,12,9,13,17,24,21,18,14) + runif(n_person, -1.5, 1.5)
			Z_f <- rnorm(n_person)
			U <- c(rmvnorm(1, sigma = add_var * phi)) 
			G1_f <- rowSums(peddrop( ped1[1:n_person, c("fa","mo")], 0.1 )) # Assumes same pedigree structure across families
			G2_f <- rowSums(peddrop( ped1[1:n_person, c("fa","mo")], maf_G2 ))
			
			Y_f <- rbinom(n_person, size = 1, prob = plogis(int + U + LG * ((G1_f > 0) & (G2_f > 0)) +	age_f * beta[1] + sex[ped1$fam == i_fam] * beta[2] + Z_f * beta[3]))
			
			if(sum(Y_f) > -1) ascert <- 1 #At least 6 affected
		}
		Y <- append(Y, Y_f)
		age <- append(age, age_f)
		Z <- append(Z, Z_f)
		G1 <- append(G1, G1_f)
	}
	fam_data = data.frame(Y=Y, age = age, sex = sex, Z = Z, idU = 1:n_total)

### MCMC glmm
	prior1 <- list(R = list(V = 1, fix=1),
	G = list(G1 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=10^2)))
	mMCMC <- MCMCglmm(Y ~ 1 + age + sex + Z, random =~idU, family = "categorical",
	ginverse = list(idU = Prec), prior=prior1,
	data = fam_data, verbose = F, nitt=250e3, thin=500, burnin=5e3)
	k <- (16*sqrt(3)/(15*pi))^2
	pvalstot <- list(betas = posterior.mode(mMCMC$Sol/sqrt(1 + k * rowSums(mMCMC$VCV))),
		vars = posterior.mode(mMCMC$VCV[,1] / (1 + k * rowSums(mMCMC$VCV))),
		vars_v2 = posterior.mode(mMCMC$VCV[,1] / (1 + k)),
		effect_size = effectiveSize(mMCMC$VCV[,1]))
	return(pvalstot)
}

M = rbind(c(.8,6.3))
n_M = dim(M)[1]
# Will be a list (over Tu/Te) of a list (over rep) of matrices (GMMAT and MCMC)
res = rep(list(list()), n_M) 
n_rep = 1000; kt <- 0 # To keep track of iterations

for(i in 1:n_M) {
  res[[i]] = mclapply(1:n_rep, runme, mc.cores = 12, Tu= M[i,1], Tot = M[i,2])
  save.image() #Save the res object
}

allfile <- list.files()
file.remove(allfile[grep("BATMAN.",allfile)])
  
body <- tail(strsplit(getwd(),"/")[[1]],1); body <- paste("Folder: ", body) 
system(paste("echo \"",body,"\" | mailx -s \"Sim done\" joelle.mbatchou@gmail.com"), ignore.stdout = T, ignore.stderr = T)