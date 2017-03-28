# 200K MCMC iterations W ASCERT 6+
library(mvtnorm)
library(parallel)
library(reshape)

library(MCMCglmm)
server <- strsplit(getwd(),"/")[[1]][2]
lib_folder <- paste0("/",server, "/mbatchou/GMMAT") 
gmmat_load <- require(GMMAT, lib.loc = lib_folder)
if (!gmmat_load) {
lib_folder <- paste0("/",server, "/mbatchou/GMMAT_compute8")
require(GMMAT, lib.loc = lib_folder)
}

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
	LG <- log(1.1) # increase of 10% in the odds of being a case
	int <- qlogis(0.1) - beta[1] * 31.86 -  beta[2] * 1.5 - LG * mean_G1G2  # Prevalence - b * E(X) - lambda * mean of epistatis indicator
	pvals_t1err = NA; pvals_power = NA
	
	while(any(is.na(pvals_t1err)) | any(is.na(pvals_power))){
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
				
				if(sum(Y_f) > 5) ascert <- 1 #At least 6 affected
			}
			Y <- append(Y, Y_f)
			age <- append(age, age_f)
			Z <- append(Z, Z_f)
			G1 <- append(G1, G1_f)
			# G2 <- append(G2, G2_f)
			# tapply(Y, (G1 > 0) & (G2 > 0), mean)
		}
		fam_data = data.frame(Y=Y, age = age, sex = sex, Z = Z, idU = 1:n_total)

### MCMC glmm
		prior1 <- list(R = list(V = 1, fix=1),
		G = list(G1 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=10^2)))
		mMCMC <- MCMCglmm(Y ~ 1 + age + sex + Z, random =~idU, family = "categorical",
		ginverse = list(idU = Prec), prior=prior1,
		data = fam_data, verbose = F, nitt=250e3, thin=500, burnin=5e3)
		pi_hat <- as.vector(predict(mMCMC, marginal = ~idU, type = "response")) # E(pi)

## Generating unassociated Gs
		n_G <- 1000 # 50 * n_rep replicates for type 1 error
		Gmat <-  cbind(G1, replicate(n_G, rowSums(peddrop( ped1[,c("fa","mo")], 0.3 )))) # G1 is first column + n * n_G matrix 
# Numerator
		T_num_G <- as.vector(crossprod(Gmat, Y - pi_hat))^2
		
# Retrospective denominator
		get_denum_retro <- function(G, Matdenum){ 
			mafs <- colMeans(G) / 2
			2 * mafs * (1 - mafs) * Matdenum
		}
		YPhiY <- crossprod(Y - pi_hat, Phimat) %*% (Y - pi_hat) #marginal mean
		T_denum_G <- get_denum_retro(Gmat, YPhiY)
		pvalstot <- pchisq(T_num_G / T_denum_G, df = 1, lower.tail = F)
		pvals_t1err <- pvalstot[-1]
		pvals_power <- pvalstot[1]

		numend <- sample(1e6,1)
		filend <- paste0(numend, '.txt')
# GMMAT
		mGMMAT <- glmmkin(Y ~ age + sex + Z,  data = fam_data, kins = Phimat, family = binomial(link="logit"))
		write.table(cbind(paste0("SNP", 1:(n_G+1)), t(Gmat)), paste0('genofile', filend), row.names=F, col.names=F, quote=F, sep='\t')
		capture.output(glmm.score.text(mGMMAT$res, mGMMAT$P, infile = paste0('genofile', filend), outfile = paste0('gmmat_out', filend)), file = "/dev/null")
		if(file.exists(paste0('gmmat_out', filend))){
		pvalstot <- as.numeric(read.table(paste0('gmmat_out', filend), header=T)[,6])
		pvals_t1err <- cbind(pvals_t1err, pvalstot[-1]) 
		pvals_power <- append(pvals_power, pvalstot[1]) 
		} else{
		pvals_t1err <- cbind(pvals_t1err, NA) 
		pvals_power <- append(pvals_power, NA) 
		write(paste('Error in GMMAT', kt), "errorfile.txt", append = T)}

# CARAT
		write.table(cbind(1,ped1[,-1],fam_data[,c(3,1,2:4)]), paste0('phenofile', filend), sep = '\t', quote=F, col.names=F, row.names=F)
		write(c('Chr','SNP','cm','bp',1:n_total), paste0('genofile', filend), ncolumns = n_total + 4, sep='\t')
		write.table(cbind(1, paste0("SNP", 1:(n_G+1)), 1,50, t(apply(Gmat,2, function(g) unlist(list(c(1,1),c(1,2),c(2,2))[g+1])))),paste0('genofile', filend), row.names=F, col.names=F,  sep='\t', quote = F, append = T)
		write.table(cbind(1,kin[kin$v>0,1:3]),paste0('grmfile', filend), row.names=F, col.names=F,sep=' ')
		com_carat <- paste('./CARAT -p',paste0('phenofile', filend),'-g',paste0('genofile', filend),'-R ',paste0('grmfile', filend),'-o',numend)
		ret <- system(com_carat, ignore.stdout = T, ignore.stderr = T)
		if ((!file.exists(paste0(numend,'_out.txt'))) | (ret !=0)){
			pvals_t1err <- cbind(pvals_t1err, NA) 
			pvals_power <- append(pvals_power, NA) 
			allfile <- list.files()
			file.remove(allfile[grep(numend,allfile)])
			write(paste('Error in CARAT', kt), "errorfile.txt", append = T)
		} else {
			if( length((pvalstot <- as.numeric(read.table(paste0(numend,'_out.txt'), header=T)[,4]))) == 0){
				pvals_t1err <- cbind(pvals_t1err, NA) 
				pvals_power <- append(pvals_power, NA) 
			} else {
				pvals_t1err <- cbind(pvals_t1err, pvalstot[-1]) 
				pvals_power <- append(pvals_power, pvalstot[1]) 
			}
		}

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
		system(paste('rm -R', numend))
		
# Erase files
		allfile <- list.files()
		file.remove(allfile[grep(numend,allfile)])
	}
	pvalstot <- list(Type1Err = pvals_t1err, Power = pvals_power)
    return(pvalstot)
}

M = 
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



########################### 
################## Analysis
plot_power <- function(power_mat){
na_i <- n_p > 0
props <- seq(0,100,by=20)[na_i]
sds <- 2*sqrt(.25/n_p[na_i])

pdf("power_plot.pdf")
errbar(props, (y<-power_mat[na_i,1]), y+sds, y-sds, type = "b", pch=16, col='red', errbar.col = 'red', ylim = c(0,1), cap = .03,
 ylab = "Power", xlab = "Proportion of variance on logit scale due to polygenic effects/covariates", xaxt = "n" )
for(i in 2:4){
	errbar(props, (y<-power_mat[na_i,i]), y+sds, y-sds, type = "b",add=T, cap = .03,pch = c(15,17,18)[i-1], col=c('blue','lightgreen','purple')[i-1], errbar.col = c('blue','lightgreen','purple')[i-1])
}
axis(1, at = props, label = c("0/100", "20/80","40/60","60/40","80/20","100/0")[na_i])
legend('topright', horiz = F, legend = c('MCMCglmm', 'GMMAT','CARAT','CERAMIC'), pch = c(16,15,17,18), col = c('red','blue','lightgreen','purple'), cex=1.4, bty="n")
dev.off()
}

pvals_analysis <- function(p_vals, al){
    n <- length(p_vals); i <- 1
    alpha <- c(.005,.01,.05)
    result <- data.frame("alpha" = alpha, "Err_rate" = NA,"SE" = NA,"p_value" = NA)
    for(a in alpha){
      result$p_value[i] <- prop.test(x = sum(p_vals < a), n, p = a)$'p.val'
      result$Err_rate[i] <- round(as.numeric(prop.test(x = sum(p_vals < a), n, p = a)$'est'),4)
      result$SE[i] <- round(sqrt(mean(p_vals < a)*(1-mean(p_vals < a)) / n),4)
      i <- i + 1
    }; result$Err_rate[alpha==al]
  }

draw_pval <- function(ind, p_vals_list, meth_list, percent){
  p_vals <- p_vals_list[[ind]]
  method <- meth_list[ind]
  n <- length(p_vals)
  uni2 <- rank(p_vals, ties.method='max'); names(uni2)=c()
  plot(-log10(uni2/(n+1)),-log10(p_vals), type="n", ylab=expression('Observed (-log'[10]*' p-value)'), xlab=expression('Expected (-log'[10]*' p-value)'))
  a=1:n
  high <- qbeta(0.025, a, rev(a))
  low <- qbeta(0.975, a, rev(a))
  polygon(-log10(c(a/n,rev(a/n))), -log10(c(high, rev(low))), col ='gray', border = NA)
  points(-log10(uni2/(n+1)),-log10(p_vals), pch=16,cex=.3) 
  title(expression(bold('Plot of observed vs. expected -log'[10]*' p-values')),line=2)
  title(paste('With', method), line=.5)
  abline(a=0,b=1,col="red")
}


################################## for 2 settings of Tu/Tot
n_p=rep(0, 6)
method_used = c('MCMCglmm', 'GMMAT','CARAT','CERAMIC')
power_mat = matrix(NA, 6, 4)
type1err <- matrix(NA, 6, 4, dimnames = list(NULL, method_used))

index <- 5:6
for(k in index){
assign(paste0('p_vals',k,'_t1err'), vector("list", 4)) # List for 4 methods
assign(paste0('p_vals',k,'_pow'), vector("list", 4)) # List for 4 methods
}

for(i in 1:4){
	for(j in 1:n_rep){
	for(k in 1:length(index)){
		if(is.double(unlist(res[[k]][[j]]))){
			eval(parse(text=paste0('p_vals',index[k],'_t1err[[',i,']] <- append(p_vals',index[k],'_t1err[[',i,']] , res[[',k,']][[',j,']]$Ty[,',i,'])')))
			eval(parse(text=paste0('p_vals',index[k],'_pow[[',i,']] <- append(p_vals',index[k],'_pow[[',i,']] , res[[',k,']][[',j,']]$Po[',i,'])')))
			if (i ==1) n_p[index[k]] = n_p[index[k]] + 1; }
			}
		}
}


index = 1:6
level = .05
# Power
for(k in 1:length(index)){
	eval(parse(text=paste0('power_mat[',index[k],',] <- unlist(lapply(p_vals',index[k],'_pow, function(x) mean(x<level)))')))
}
library(Hmisc)
plot_power(power_mat)
# Type 1 Error
props <- seq(0,100,by=20)
for(k in 1:length(index)){
	eval(parse(text=paste0('type1err[',index[k],',] <- unlist(lapply(p_vals',index[k],'_t1err, pvals_analysis, al = .05))')))
}
type1err; .05+c(-1,1)*2*sqrt(.05*(1-.05)/5e4)
png("type_1_error_plots%01d.png", width = 12, height = 9, units = 'in', res=200)
par(mfrow = c(2,2))
for(k in 1:length(index)){
	eval(parse(text=paste0('lapply(1:4, draw_pval, p_vals',index[k],'_t1err, method_used)')))
	title(paste('Polygenic effect / Covariates:',props[index[k]],'/',100-props[index[k]]), outer=TRUE, cex=1.1, line = -1.5)
}
dev.off()


################################## for 5 settings of Tu/Tot
n_p=rep(0, n_M)
method_used = c('MCMCglmm', 'GMMAT','CARAT','CERAMIC')
power_mat <- matrix(NA, n_M, 4)
type1err <- matrix(NA, n_M, 4, dimnames = list(NULL, method_used))

for(k in 1:n_M){
assign(paste0('p_vals',k,'_t1err'), vector("list", 4)) # List for 4 methods
assign(paste0('p_vals',k,'_pow'), vector("list", 4)) # List for 4 methods
}

for(i in 1:4){
	for(j in 1:n_rep){
	for(k in 1:n_M){
		if(is.double(unlist(res[[k]][[j]]))){
			eval(parse(text=paste0('p_vals',k,'_t1err[[',i,']] <- append(p_vals',k,'_t1err[[',i,']] , res[[',k,']][[',j,']]$Ty[,',i,'])')))
			eval(parse(text=paste0('p_vals',k,'_pow[[',i,']] <- append(p_vals',k,'_pow[[',i,']] , res[[',k,']][[',j,']]$Po[',i,'])')))
			if (i ==1) n_p[k] = n_p[k] + 1; }
			}
		}
}

level <- .05
# Power
for(k in 1:n_M){
	eval(parse(text=paste0('power_mat[',k,',] <- unlist(lapply(p_vals',k,'_pow, function(x) mean(x<level)))')))
}
library(Hmisc)
plot_power(power_mat)

# Type 1 Error
props <- seq(0,100,by=20)
for(k in 1:n_M){
	eval(parse(text=paste0('type1err[',k,',] <- unlist(lapply(p_vals',k,'_t1err, pvals_analysis, al = .05))')))
}
type1err; 2*sqrt(.05*(1-.05)/n_p)
png("type_1_error_plots%d.png", width = 12, height = 9, units = 'in', res=200)
par(mfrow = c(2,2))
for(k in 1:n_M){
	eval(parse(text=paste0('lapply(1:4, draw_pval, p_vals',k,'_t1err, method_used)')))
	title(paste('Polygenic effect / Covariates:',props[k],'/',100-props[k]), outer=TRUE, cex=1.1, line = -1.5)
}
dev.off()
