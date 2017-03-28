# Quantitative Yperm
# With VT
# No ascertainment
# High prev setting (30-33%)
# Vary add_var (0-20% var(pi)/var(Y))

pckg_list <- c("mvtnorm","MCMCglmm", "parallel", "matrixStats", "Matrix", "data.table","reshape")
dummy <- lapply(pckg_list, require, character.only = TRUE)
# lib = "../R_lib_compute78/")

# Parameters for simulation
N_fam <- 45
n_person <- 22
n_asc <- 0 # Degree of ascertainment
n_total <- N_fam * n_person
method_name <- "VT"   # VT or ws
Yperm_type <- "quant"  #quant or bin
n_perm <- 20e3
nsets <- 150 # Number of RV sets/realization in analysis
n_marker <- 50 # Number of rvs/set
ncov <- 4

# Store haplo info
haplo_file <- fread("../haplo.txt")

# Make pedigree
ped1 <- data.table(fam = NA, ind = 1:n_total, fa = NA, mo = NA)
fa_ind <- c(NA,NA,NA,1,1,NA,1,NA,1,NA,4,4,4,6,6,6,8,8,10,10,10,10)
mo_ind <- c(NA,NA,NA,2,2,NA,2,NA,2,NA,3,3,3,5,5,5,7,7,9,9,9,9)
ped1[,`:=`(mo = (rep(mo_ind, N_fam) + n_person*rep(0:(N_fam-1), each = n_person)),
fa = (rep(fa_ind, N_fam)+ n_person*rep(0:(N_fam-1), each = n_person)), fam = rep(1:N_fam, each = n_person))]
ped1[is.na(ped1)] <- 0
write.table(ped1,'pedigree',row.names=FALSE, col.names=FALSE,sep='\t')
sex <- rep(c(1,2,2,1,2,1,2,1,2,1,2,1,1,2,1,2,1,2,2,1,2,1),N_fam) #constant

# Generate relatedness data
phi <- as.matrix(fread("add_phi_mat_n22")); colnames(phi) <- 1:n_person
setnames(kin <- data.table(melt(diag(N_fam)%x%phi)), 1:3, c("person2","person1","kincoef"))
# setcolorder(kin, c("person1","person2","kin"))
kin <- kin[person1<=person2 & ((person1-1)%/%n_person) == ((person2-1)%/%n_person)]
kin[, fam := (person1-1)%/%n_person+ 1]
# For CERAMIC -- same for all realizations
kin[kincoef == 1, kincoef:= 0]; kin[, kincoef := kincoef/2]
write.table(kin[,.(fam, person1, person2, kincoef)],'grmfile', row.names=F, col.names=F,sep=' ')

dummy <- suppressMessages(crossprod( Matrix(diag(2) %x% phi), diag(2*n_person)))
dummy <- suppressMessages(crossprod(tcrossprod(Diagonal(2) %x% chol(phi), diag(sqrt(1:(2*n_person)))), eigen(dummy)$ve)); rm(dummy) ## Suppress all future warnings msgs from Matrix pckg

get_vt_stat <- function(ped1, pheno, mu_y, delt, GCV_mat){
	if(Yperm_type == "bin"){
		ncaseobs <- sum(ncaseobsF <- tapply(pheno, ped1[,fam], sum))
		ytot <- cbind(pheno, sapply(1:n_perm, function(i){
			repeat{
				Yp <- mu_y + GCV_mat %*% sample(delt)
				thresh <- quantile(Yp, 1 - ncaseobs/n_total)
				Yp <- as.numeric(Yp > thresh)
				dev <- abs(tapply(Yp, ped1[,fam], sum) - ncaseobsF)
				# Set overall maximum number of deviations to 10%
				# Set maximum allowable deviation in ncase within families to 20% (4.4 per 22 inds)
				if((sum(dev) < .1 * n_total) & all(dev< .2 * n_person) & (sum(Yp) == ncaseobs)) break;
			}
			return(Yp)
		}))
	} else ytot <- cbind(pheno, mu_y + tcrossprod(GCV_mat, t(replicate(n_perm, sample(delt)))))
	ytot <- scale(ytot, center = TRUE, scale = FALSE)
	varsY <- apply(ytot, 2, var)
	tol <- 1e-6
	
	sapply(seq(nsets), function(i_set) {
		genos <- genosim(ped1[,.(fa,mo)], n_marker)
		mafs <- (1 + colSums(genos, na.rm = TRUE))/(2 + 2 * nrow(genos))
		t.maf <- c(sort(unique(mafs)),1)
		groups <- as.numeric(cut(mafs, t.maf, right = F))
		z <- lapply(seq(max(groups)), function(i_T){
			genos_sub <- genos[,groups == i_T]
			denum <- sum(genos_sub^2)
			num <-  rowSums(crossprod(ytot, genos_sub))
			list(num = num, denum = denum)
		})
		nums <- t(apply(sapply(z, function(a) a$num), 1, cumsum))
		denums <- outer(varsY, cumsum(sapply(z, function(a) a$denum)))
		z <- rowMaxs(abs( nums/sqrt(denums) ))
		eq <- sum( abs(z[-1]-z[1]) < tol )
		num_pval <- sum((z[-1] - z[1]) > tol)+ sum(runif(eq) > .5)
		(num_pval + 1.0)/ length(z)
	})
}

get_ws_stat <- function(ped1, pheno, mu_y, delt, GCV_mat){
	tol <- 1e-6
	ncaseobs <- sum(ncaseobsF <- tapply(pheno, ped1[,fam], sum))
	# convert Yp to binary
	ytot <- cbind(pheno, sapply(1:n_perm, function(i){
		accept = 0
		while(accept == 0){
			Yp <- mu_y + GCV_mat %*% sample(delt)
			thresh <- quantile(Yp, 1 - ncaseobs/n_total)
			Yp <- as.numeric(Yp > thresh)
			dev <- abs(tapply(Yp, ped1[,fam], sum) - ncaseobsF)
			# Set overall maximum number of deviations
			# Set maximum allowable deviation in ncase within families
			if((sum(dev) < .05 * n_total) & all(dev< .3 * n_person) & (sum(Yp) == ncaseobs)) accept = 1
		}
		return(Yp)
	}))
	indy <- ytot==0
	
	sapply(seq(nsets), function(a) {
		genos <- genosim(ped1[,.(fa,mo)], n_marker)
		w <- sapply(1:ncol(indy), function(i){
			y <- indy[,i]  # 1 if cont and 0 if case
			w <- (1 + colSums(genos[y,]))/(2 + 2 * (n_total-ncaseobs))
		})
		w <- 1 / sqrt(n_total * w * (1-w))
		Rscore <- apply(genos %*% w, 2, rank) # score for each person based on Y replicate
		z <- colSums((!indy)*Rscore)
		
		eq <- sum( abs(z[-1]-z[1]) < tol )
		num_pval <- sum((z[-1] - z[1]) > tol)+ sum(runif(eq) > .5)
		(num_pval + 1.0)/ ncol(indy)
	})
}

gendrop <- function (parentmat, MAF) { # parentmat is matrix of parental info
	nfound <- parentmat[fa == 0 & mo == 0, .N]
	npeople <- nrow(parentmat)
	haplos <- matrix(0, npeople, 2) # Assign haplotype to each individual in the pedigree
	haplos[(parentmat[, fa] == 0) & (parentmat[, mo] == 0),] <- matrix(1:(2 * nfound), nfound, 2, byrow=T)

	for (i in 1:npeople) {
		#/* When both parents are in the pedigree drop alleles down the pedigree */
		if (sum(parentmat[i]) != 0) {
			#/* Randomly choose chrom. 0 or 1 to transmit to kid at each locus
			trchrom = sample(1:2, 2, replace = T)
			#// For GT haplotype
			haplos[i, 1] <- haplos[parentmat[i, fa], trchrom[1]];
			haplos[i, 2] <- haplos[parentmat[i, mo], trchrom[2]];
		} 
	}

	state <- as.numeric(runif(2 * nfound) < MAF) #1 if minor, 0 if major allele
	genos <- state[haplos]
	genos <- genos[1:npeople] + genos[-(1:npeople)] 
	return(genos) # Vector of MAC counts
}

genosim <- function (parentmat, n_marker) { # Returns GT info for one rv set of n_markers
	nfound <- sum(is_found <- parentmat[,fa + mo] == 0)
	npeople <- nrow(parentmat)
	haplos <- matrix(0, npeople, 2) # Assign haplotype to each individual in the pedigree
	haplos[is_found,] <- matrix(1:(2 * nfound), nfound, 2, byrow=T)

	trchrom = matrix(sample(1:2, 2 * (npeople - nfound) , replace = T), (npeople - nfound), 2)
	for (i in (1:npeople)[!is_found]) {
		#/* When both parents are in the pedigree drop alleles down the pedigree */
			#/* Randomly choose chrom. 0 or 1 to transmit to kid at each locus
			haplos[i,] <- sapply(1:2, function(a) haplos[parentmat[i][[a]], trchrom[which((1:npeople)[!is_found] == i), a]])
	}

	# Get founder haplotypes 
	founder_haplo <- sample(seq(10000), 2 * nfound, replace = F)
	state <- -(haplo_file[founder_haplo]-2) #1 if minor, 0 if major allele
	genos <- state[c(haplos)] 
	genos <- genos[1:npeople] + genos[-(1:npeople)] 
	genos <- genos[,which(unlist(lapply(genos, function(x) sd(x) > 0))), with=F]
	
	if(ncol(genos) < n_marker) return(NA)
	# Scramble which GT column is picked for first n_marker polymorphic genes
	GTpicked <- sample(1:ncol(genos))[1:n_marker]
	return(as.matrix(genos[,GTpicked, with = F])) # n_people * n_marker matrix of rvs
}

runme <- function (dummy, Tu, Tot)
{
	write(dummy, "iteration_track.txt")
	add_var <- Tu * Tot
	# betas <- c(.05, 1, .5)
	# old_v_xb <- sum(c(363.6669,.25,1) * betas^2)
	# S_XB <- sqrt(Tot * (1 - Tu) / old_v_xb) # Ratio of target variance to var(X)
	# betas <- betas * S_XB
	S_XB <- sqrt(Tot * (1 - Tu) / (363.67+.25+1)) # Wrong way of fixinf Var(Xb) -- used to compare to Param boot results
	betas <- c(.05, 1, .5) * S_XB
	LG <- log(1) # NONE -- increase of 10% in the odds of being a case
	
	# Covariates in generating model
	int <- qlogis(0.3) - sum(betas * c(31.86, 1.5, 0)) - LG * .1425  # Prevalence - b * E(X) - lambda * .1425
	# G1 <- gendrop( ped1[, .(fa,mo)], 0.1 ) # Assumes same pedigree structure across families
	# G2 <- gendrop( ped1[, .(fa,mo)], 0.5 )
	U <- c(t(rmvnorm(N_fam, sigma = add_var * phi)))
	age <- rep(c(73,75,46,43,40,46,40,43,47,51,18,21,15,15,12,9,13,17,24,21,18,14), N_fam) + runif(n_total, -1.5, 1.5)
	Z <- rnorm(n_total)
	###--###
	covars <- int + U  +	cbind(age, sex, Z) %*% betas #+ LG * ((G1 > 0) & (G2 > 0)) 

	Y <- c(sapply(1:N_fam, function(i_fam){
		Y_f <- -5
		while(sum(Y_f) < n_asc) Y_f <- rbinom(n_person, size = 1, prob = plogis(covars[1:n_person + (i_fam - 1) * n_person]))
		return(Y_f)
	}))
	
	fam_data <- data.table(ped1, sex, Y, age, Z)
	numend <- paste0(dummy,"_",sample(1e6,1))
	filend <- paste0(numend, '.txt')
	
	# CERAMIC
	system(paste('mkdir', numend))	
	write.table( fam_data[,.(fam, ind, Y+1, age, sex, Z)], paste0('./',numend,'/','phenofile'), sep = '\t', quote=F, col.names=F, row.names=F)
	com_ceramic <- paste0('(cd ./',numend,' && exec ../BATMAN -p phenofile -k ../grmfile -N)')
	ret <- system(com_ceramic, ignore.stdout = T, ignore.stderr = T)
	if ((!file.exists(paste0('./',numend,'/BATMANtest.phenoestimates'))) | (ret !=0)){
		write('Error in CERAMIC', "errorfile.txt", append = T)  
		system(paste0('rm -R *', numend, "*"))
		return(result_out <- list(Type1Error = -1, herit = NA))
	}
	pvalstot <- fread(paste0('./',numend,'/BATMANtest.phenoestimates'), skip = 5, header=T)
	xi <- pvalstot[1, Estimate]
	est_beta <- pvalstot[4:7, Estimate]
	system(paste0('rm -R *', numend, "*"))

	## mu
	mu_y <- drop(plogis((X <- model.matrix(~age+sex+Z, fam_data)) %*% est_beta))
	## gamma
	gam <- mu_y * (1-mu_y)
	## sigma_f
	sigf <- xi * phi + (1-xi) * diag(n_person)
	C_sig <- Diagonal(N_fam) %x% chol(sigf)
	## w
	w <- backsolve(C_sig, crossprod(diag(sqrt(gam)), X), transpose=TRUE)
	## Omeg
	suppressMessages(Omeg <- Diagonal(n_total) - w %*% solve(crossprod(w), t(w)))
	## V
	omeg_eig <- eigen(Omeg, symmetric = T); V <- omeg_eig$ve[,omeg_eig$va>.9]
	## VCe
	yres <- backsolve(C_sig, crossprod(diag(1/sqrt(gam)), fam_data[,Y] - mu_y), transpose=TRUE)
	delt <- drop(crossprod(V, yres))
	## gam*C*V
	suppressMessages({pre_mult_mat <- crossprod(tcrossprod(C_sig, diag(sqrt(gam))) ,V)})
	# # Check estimating eqns
	# sig_inv <- diag(N_fam) %x% solve(sigf)
	# t(X) %*% diag(sqrt(mu_y * (1-mu_y))) %*% sig_inv %*% diag(1/sqrt(mu_y * (1-mu_y))) %*% (fam_data[,Y] - mu_y)
	# mean(pre_mult_mat %*% sample(delt))
	rm(X,gam,sigf,C_sig,w,Omeg,omeg_eig, V, yres)
	
	# Simulate GT, compute VT & get pvals
	if(method_name == "VT"){
		pvalsEE <- get_vt_stat(ped1, fam_data[,Y], mu_y, delt, pre_mult_mat)
	} else pvalsEE <- get_ws_stat(ped1, fam_data[,Y], mu_y, delt, pre_mult_mat)
	
	return(list(Type1Error = pvalsEE, herit = xi))
}

############################ 
################### Analysis
# M = data.frame(Tu = seq(0,1,.2), Tot = c(1.1,1.2,1.2,1.3,1.3,1.3))
M = data.frame(Tu = 0, Tot = 0)
n_M = nrow(M)
# Will be a list (over Tu/Te) of a list (over rep) of matrices (GMMAT and MCMC)
res = rep(list(list()), n_M) 
n_rep = 150

for(i in 1:n_M) {
	res[[i]] = mclapply(1:n_rep, runme, mc.cores = 16, Tu= M$Tu[i], Tot = M$Tot[i])
	save.image() #Save the res object
}

pvals <-lapply(res, function(x) unlist(lapply(x, function(a) if(length(a) > 1) a$Ty)))
herit <-lapply(res, function(x) unlist(lapply(x, function(a) if(length(a) > 1) a$her)))

get_type1errs <- function(p){
	# summary(p)
	n=length(p); print(paste("n =",n))
	
	# hist(p, main="Histogram of permutation-based p-values",freq = F,breaks=100)
	# abline(h=1,col="red")

	result <- data.frame("alpha"= c(.005,.01,.05), "Err_rate"= NA, "SE"= NA, "p_value"= NA)
	result <- t(sapply(1:3, function(i){
		a <- result$alpha[i]
		result$p_value[i] <- round(prop.test(x = sum(p < a), n, p = a)$'p.val',4)
		result$SE[i]=round(sqrt(mean(p<a)*(1-mean(p<a))/n),4)
		result$Err_rate[i]=round(as.numeric(prop.test(x = sum(p<a),n,p = a)$'est'),4)
		result[i,]
	}))
	print(result)
	result
}

x <- unlist(get_type1errs(pvals[[1]])[,4])

body <- tail(strsplit(getwd(),"/")[[1]],1); body <- paste("Folder: ", body , "\n", "pvalue at .05 = ", x[3],"\n", "pvalue at .01 = ", x[2],"\n", "pvalue at .005 = ", x[1])
system(paste("echo \"",body,"\" | mailx -s \"Sim done\" joelle.mbatchou@gmail.com"), ignore.stdout = T, ignore.stderr = T)



