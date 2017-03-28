# Changes:
# Ascertainment present (low preval.)
# Decreased by 1/2 number of MCMC iterations
# Exact computation of polygen. add. var. adjusting for iid residuals in MCMCglmm

library(mvtnorm)
library(MCMCglmm)

# Integrate out the iid residuals
normal.logistic <- function(mu, v) {
	int.foo <- function(x, mu, v) {
		plogis(x) * dnorm(x, mu, sqrt(v))
	}
	integrate(int.foo, qnorm(0.0001, mu, sqrt(v)), 
	qnorm(0.9999, mu, sqrt(v)), mu, v)[[1]]
}

runme = function (dummy, famsize)
{
	N_fam <- famsize;
	n_person <- 22 ;
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

# Generate data
	phi <- as.matrix(read.table("add_phi_mat_n22")); colnames(phi) <- 1:n_person
	phi_inv <- solve(phi) #22x22 matrix
	sex <- c(1,2,2,1,2,1,2,1,2,1,2,1,1,2,1,2,1,2,2,1,2,1)

	add_v <- 5.8 * .2
	beta<- sqrt(5.8 * (1 - .2) / .25) 
	int <- qlogis(0.025) - beta * 1.5
	for(i_fam in 1:N_fam){
		while(sum(Y_f <- 0) < 4){
			U <- rmvnorm(1, sigma = add_v * phi)
			Y_f <- rbinom(n_person, size = 1, prob = plogis(int + beta * sex + U))
		}
		Y <- append(Y, Y_f)
	}
	fam_data = data.frame(Y= Y, X= rep(sex, N_fam) , idu = 1:n_total)

##MCMC glmm
	prior1 <- list(R = list(V = 1, fix=1),
	G = list(G1 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=10^2)))
	mMCMC<-MCMCglmm(Y ~ 1 + X, random =~idu, family = "categorical",
	ginverse = list(idu = Prec), prior=prior1,
	data = fam_data, verbose = F, nitt=150e3, thin=500, burnin=15e3)
# ck <- ((16*sqrt(3))/(15*pi))^2
	m <- sqrt(c(posterior.mode(mMCMC$VCV[,1]), colMeans(mMCMC$VCV)[1]))
	m_est <- qlogis(sapply(m, normal.logistic, v= 1))^2
	effS <- effectiveSize(mMCMC$VC[,1])
	return (c(m_est, effS))
}


M = c(50,100,150)
n_M = length(M)
res = rep(list(list()), n_M)
n = 150
target = 5.8*.2

library(parallel)
for(i_famsiz in 1:n_M) {
	res[[i_famsiz]] = mclapply(1:n, runme, mc.cores = 15, famsize=M[i_famsiz])
}

body <- tail(strsplit(getwd(),"/")[[1]],1); body <- paste("Folder: ", body) 
system(paste("echo \"",body,"\" | mailx -s \"Sim done\" joelle.mbatchou@gmail.com"), ignore.stdout = T, ignore.stderr = T)

# pdf("mode_mean_est_var_MCMC_reduced.pdf", height=4, width=5.5)
# est_vec <- rep(list(matrix(NA, n, 2)), n_M); names(est_vec) = paste0("n_", M)
# for(i in 1:n_M) {
# for(j in 1:n)  {
		# est_vec[[i]][j,] <- res[[i]][[j]]
# }
# }

# sapply(seq(n_M), function(i){ boxplot(est_vec[[i]], col = rep(c(3,5)),
	# main = paste("Posterior estimates for the variance for MCMCglmm, nf=",M[i]), ylim=c(0,8/i))
# legend('topright', legend = c("Mode", "Mean"), col=c(3,5), bty = "n", pch =16, cex=.8)
# abline(h = target, lwd=2, lty = 2, col=4)
# })
# dev.off()

# lapply(seq(n_M), function(i) apply(est_vec[[i]], 2, function(v) mean(v)+c(-1,1)*qnorm(1-.05/2)*sd(v)/sqrt(length(v)) ) )

