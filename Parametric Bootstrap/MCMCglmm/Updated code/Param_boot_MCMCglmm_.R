# Changes:

pckg_list <- c("mvtnorm","MCMCglmm", "parallel", "matrixStats", "Matrix", "data.table","reshape")

server <- strsplit(getwd(),"/")[[1]][2]
gmmat_load <- tryCatch(suppressWarnings(library(GMMAT)), error = function(dummy){
	lib_folder <- paste0("/",server, "/mbatchou/GMMAT_compute78")
	library(GMMAT, lib.loc = lib_folder)
	return(0)}
)

dyn.load("ParamBoot")

# Parameters for simulation
N_fam <- 45
n_person <- 22
n_asc <- 0 # Degree of ascertainment
n_total <- N_fam * n_person
method_name <- "VT"   # VT or ws
n_perm <- 20e3
nsets <- 150 # Number of RV sets/realization in analysis
n_marker <- 50 # Number of rvs/set
ncov <- 4

# Make pedigree
ped1 <- data.table(fam = NA, ind = 1:n_total, fa = NA, mo = NA)
fa_ind <- c(NA,NA,NA,1,1,NA,1,NA,1,NA,4,4,4,6,6,6,8,8,10,10,10,10)
mo_ind <- c(NA,NA,NA,2,2,NA,2,NA,2,NA,3,3,3,5,5,5,7,7,9,9,9,9)
ped1[,`:=`(mo = (rep(mo_ind, N_fam) + n_person*rep(0:(N_fam-1), each = n_person)),
fa = (rep(fa_ind, N_fam)+ n_person*rep(0:(N_fam-1), each = n_person)), fam = rep(1:N_fam, each = n_person))]
Prec <- inverseA(ped1[,c('ind','mo','fa')])$Ainv
ped1[is.na(ped1)] <- 0
write.table(ped1,'pedigree',row.names=FALSE, col.names=FALSE,sep='\t')
sex <- rep(c(1,2,2,1,2,1,2,1,2,1,2,1,1,2,1,2,1,2,2,1,2,1),N_fam) #constant


# Generate data
phi <- as.matrix(read.table("add_phi_mat_n22")); colnames(phi) <- 1:n_person
Phimat <- diag(N_fam)%x%phi
setnames(kin <- data.table(melt(Phimat)), 1:3, c("person2","person1","kincoef"))
# setcolorder(kin, c("person1","person2","kin"))
kin <- kin[person1<=person2 & ((person1-1)%/%n_person) == ((person2-1)%/%n_person)]
kin[, fam := (person1-1)%/%n_person+ 1]
write.table(kin[,.(fam, person1, person2, kincoef)],'kinship', row.names=F, col.names=F,sep=' ')
sex <- rep(c(1,2,2,1,2,1,2,1,2,1,2,1,1,2,1,2,1,2,2,1,2,1),N_fam) #constant

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
	LG <- log(1) # increase of 10% in the odds of being a case

	int <- qlogis(0.3) - sum(betas * c(31.86, 1.5, 0)) - LG * 2 * .1  # Prevalence - b * E(X) - lambda * mean of epistatis indicator
	U <- c(t(rmvnorm(N_fam, sigma = add_var * phi)))
	age <- rep(c(73,75,46,43,40,46,40,43,47,51,18,21,15,15,12,9,13,17,24,21,18,14), N_fam) + runif(n_total, -1.5, 1.5)
	Z <- rnorm(n_total)
	G1 <- peddrop(ped1[, .(fa,mo)], 0.1) # 
	###--###
	covars <- int + U  +	cbind(age, sex, Z) %*% betas + LG * G1
	
	pvalstot <- NULL
	while(is.null(pvalstot)){
		Y <- c(sapply(1:N_fam, function(i_fam){
			Y_f <- -5
			while(sum(Y_f) < n_asc) Y_f <- rbinom(n_person, size = 1, prob = plogis(covars[1:n_person + (i_fam - 1) * n_person]))
			return(Y_f)
		}))
		
		fam_data <- data.table(ped1, sex, Y, age, Z, idU = 1:n_total)
### MCMC glmm
		prior1 <- list(R = list(V = 1, fix=1),
		G = list(G1 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=10^2)))
		mMCMC <- MCMCglmm(Y ~ 1 + age + sex + Z, random =~idU, family = "categorical",
		ginverse = list(idU = Prec), prior=prior1,
		data = fam_data, verbose = F, nitt=250e3, thin=500, burnin=10e3, pr = TRUE)
		ck <- ((16*sqrt(3))/(15*pi))^2
		effS <- effectiveSize(mMCMC$VCV)[1]
		est_beta <- colMeans(as.matrix(mMCMC$Sol[,1:ncov]) /sqrt(1 + ck))#niter*ncov
		est_u <- as.matrix(mMCMC$Sol[,-(1:ncov)]) /sqrt(1 + ck)
		est_s2 <- drop(mMCMC$VC[,1])/(1 + ck)
		
		# With posterior mean
		est_add_var <- mean(est_s2)
		pvalstot <- .Call("VT_analysis", Y = as.double(Y), age = as.double(age),
		sex = as.double(sex), Z = as.double(Z),
		add_var_hat = as.double(est_add_var), beta_hat = as.double(est_beta))
		# With posterior mode
		est_add_var <- posterior.mode(est_s2)
		pvalstot2 <- .Call("VT_analysis", Y = as.double(Y), age = as.double(age),
		sex = as.double(sex), Z = as.double(Z),
		add_var_hat = as.double(est_add_var), beta_hat = as.double(est_beta))
	}
	result_out <- list(pvals = cbind(pvalstot, pvalstot2), a_var_est = est_add_var)

	return(result_out)
}

############################ 
################### Analysis
M = rbind(c(.4,6.2))
n_M = dim(M)[1]
# Will be a list (over Tu/Te) of a list (over rep) of matrices (GMMAT and MCMC)
res = rep(list(list()), n_M) 
n_rep = 500; kt <- 0 # To keep track of iterations

for(i in 1:n_M) {
res[[i]] = mclapply(1:n_rep, runme, mc.cores = 17, Tu= M[i,1], Tot = M[i,2])
save.image() #Save the res object
}

body <- tail(strsplit(getwd(),"/")[[1]],1); body <- paste("Folder: ", body) 
system(paste("echo \"",body,"\" | mailx -s \"Sim done\" joelle.mbatchou@gmail.com"), ignore.stdout = T, ignore.stderr = T)


# ### analysis
# p_valsMCMC <- vector("list", n_M) # List of size n_M
# est_varMCMC <- vector("list", n_M)

# for(i in 1:n_M) {
# for (j in 1:n_rep){
	# if (! any(res[[i]][[j]][[1]] < 0)) p_valsMCMC[[i]] <- append(p_valsMCMC[[i]], res[[i]][[j]][[1]])
	# est_varMCMC[[i]] <- append(est_varMCMC[[i]], res[[i]][[j]][[2]])
# }
# }

# pdf("pvals_1pct.pdf")

# for(i in 1:n_M) {
# print(paste("Proportion of add_var : ", M[i,1]*100, "%", sep = ""))
# summary(p_valsMCMC[[i]])
# hist(p_valsMCMC[[i]], main="Using MCMCglmm estimates",freq = F, breaks=100, xlab="p-value")
# title(bquote(sigma[a]^2 ~ "/ [Var(X"*beta*") +"~ sigma[a]^2~"] ="~.(M[i,1]*100)*"%"),line=.7, cex.main = .9)
# title(bquote("n ="~.(length(p_valsMCMC[[i]]))),line=-1, cex.main = .7)
# abline(h=1,col="red")
# }

# boxplot(est_varMCMC, xaxt="n", ylab="Additive variance estimates",
		# xlab=bquote(sigma[a]^2 ~ "/ [Var(X"*beta*") +"~ sigma[a]^2~"]   (in %)"), 
		# main = paste("Variance estimates over",n_rep,"Y realizations"),
		# col = rep(c("turquoise","red"), each = n_M))
# axis(1, at=1:n_M, labels=M[,1]*100)
# abline(h= M[,1]*M[,2], col= rep("black",n_M), lty =2)
# dev.off()

# pvals_analysis <- function(p_vals){
# n <- length(p_vals); i <- 1
# alpha <- c(.005,.01,.05)
# result <- data.frame("alpha" = alpha, "Err_rate" = NA,"SE" = NA,"p_value" = NA)
# for(a in alpha){
	# result$p_value[i] <- prop.test(x = sum(p_vals < a), n, p = a)$'p.val'
	# result$Err_rate[i] <- round(as.numeric(prop.test(x = sum(p_vals < a), n, p = a)$'est'),4)
	# result$SE[i] <- round(sqrt(mean(p_vals < a)*(1-mean(p_vals < a)) / n),4)
	# i <- i + 1
# }; result
# }
# lapply(p_valsMCMC, pvals_analysis)

