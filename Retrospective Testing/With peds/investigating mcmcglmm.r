# No ascert
Y <- NULL; age <- NULL; Z <- NULL; G1 <- NULL  
		# G2 <- NULL  
		# Generate family by family
		for (i_fam in 1:N_fam){
			ascert = 0
			while(ascert == 0){
				age_f <- c(73,75,46,43,40,46,40,43,47,51,18,21,15,15,12,9,13,17,24,21,18,14) + runif(n_person, -1.5, 1.5)
				# Z_f <- rnorm(n_person)
				U <- c(rmvnorm(1, sigma = add_var * phi)) 
				
				Y_f <- rbinom(n_person, size = 1, prob = plogis(int + U + 	age_f * beta[1] + sex[ped1$fam == i_fam] * beta[2]))
				
				if(sum(Y_f) > -1) ascert <- 1 #At least 6 affected
			}
			Y <- append(Y, Y_f)
			age <- append(age, age_f)
		}
		fam_data = data.frame(Y=Y, age = age, sex = sex, idU = 1:n_total)

		
## No relatedness & no asc.

Y <- NULL; age <- NULL; Z <- NULL; G1 <- NULL  
		# G2 <- NULL  
		# Generate family by family
		for (i_fam in 1:N_fam){
			ascert = 0
			while(ascert == 0){
				age_f <- c(73,75,46,43,40,46,40,43,47,51,18,21,15,15,12,9,13,17,24,21,18,14) + runif(n_person, -1.5, 1.5)
				U <- c(rmvnorm(1, sigma = add_var * diag(n_person))) 
				
				Y_f <- rbinom(n_person, size = 1, prob = plogis(int + U + 	age_f * beta[1] + sex[ped1$fam == i_fam] * beta[2]))
				
				if(sum(Y_f) > -1) ascert <- 1 #At least 6 affected
			}
			Y <- append(Y, Y_f)
			age <- append(age, age_f)
		}
		fam_data2 = data.frame(Y=Y, age = age, sex = sex, idU = 1:n_total)

		
## No relatedness & W asc.
Y <- NULL; age <- NULL; Z <- NULL; G1 <- NULL  
		# G2 <- NULL  
		# Generate family by family
		for (i_fam in 1:N_fam){
			ascert = 0
			while(ascert == 0){
				age_f <- c(73,75,46,43,40,46,40,43,47,51,18,21,15,15,12,9,13,17,24,21,18,14) + runif(n_person, -1.5, 1.5)
				U <- c(rmvnorm(1, sigma = add_var * diag(n_person))) 
				
				Y_f <- rbinom(n_person, size = 1, prob = plogis(int + U + 	age_f * beta[1] + sex[ped1$fam == i_fam] * beta[2]))
				
				if(sum(Y_f) > 5) ascert <- 1 #At least 6 affected
			}
			Y <- append(Y, Y_f)
			age <- append(age, age_f)
		}
		fam_data3 = data.frame(Y=Y, age = age, sex = sex, idU = 1:n_total)

## Relatedness & W asc.
Y <- NULL; age <- NULL; Z <- NULL; G1 <- NULL  
		# G2 <- NULL  
		# Generate family by family
		for (i_fam in 1:N_fam){
			ascert = 0
			while(ascert == 0){
				age_f <- c(73,75,46,43,40,46,40,43,47,51,18,21,15,15,12,9,13,17,24,21,18,14) + runif(n_person, -1.5, 1.5)
				U <- c(rmvnorm(1, sigma = add_var * phi)) 
				
				Y_f <- rbinom(n_person, size = 1, prob = plogis(int + U + 	age_f * beta[1] + sex[ped1$fam == i_fam] * beta[2]))
				
				if(sum(Y_f) > 5) ascert <- 1 #At least 6 affected
			}
			Y <- append(Y, Y_f)
			age <- append(age, age_f)
		}
		fam_data4 = data.frame(Y=Y, age = age, sex = sex, idU = 1:n_total)










###############
### MCMC glmm
prior1 <- list(R = list(V = 1, fix=1),
G = list(G1 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=10^2)))
mMCMC <- MCMCglmm(Y ~ 1 + age + sex , random =~idU, family = "categorical",
ginverse = list(idU = Prec), prior=prior1,
data = fam_data, verbose = F, nitt=250e3, thin=500, burnin=5e3)


### MCMC glmm
ck <- ((16*sqrt(3))/(15*pi))^2
posterior.mode(mMCMC_noasc$Sol/sqrt(1 + ck * rowSums(mMCMC$VCV[,])))
posterior.mode(mMCMC$VCV[,1]/(1 + ck * mMCMC$VCV[,2]))


beta_ests <- matrix(NA, n_rep, length(res[[1]][[1]]$betas) )
add_var_est <- rep(NA, n_rep)
eff_siz <- rep(NA, n_rep)
for(i in 1:n_rep){
	beta_ests[i,] <- res[[1]][[i]]$bet
	add_var_est[i] <- res[[1]][[i]]$var
	eff_siz[i] <- res[[1]][[i]]$eff
}

add_var <- prod(M)
betas <- c(.05, 1, .5)
old_v_xb <- sum(c(363.6669,.25,1)*betas^2)
S_XB <- sqrt(M[2] * (1 - M[1]) / old_v_xb) # Ratio of target variance to var(X)
betas <- betas * S_XB

# # add_var_est2 is the correct way of computing add_var est.
# colMeans(beta_ests[,-1]); beta
# boxplot(list(add_var_est1,add_var_est2), ylim=c(0,20), main = "VC estimates")
# abline(h=add_var, col='red')

boxplot(eff_siz, title = "MCMC effective sample size for variance parameter")

lim_e <- 50
cond1<- c("small" ,"large") [(eff_siz>lim_e)+1];table(cond1)
boxplot(add_var_est~cond1, ylim=c(0,20), main = paste("VC estimates (eff. size >",lim_e), xlab = "Effective size")
abline(h=add_var, col='red')


