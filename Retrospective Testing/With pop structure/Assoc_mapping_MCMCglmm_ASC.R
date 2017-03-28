library(parallel)
library(reshape)
library(matrixStats)

library(MCMCglmm)
server <- strsplit(getwd(),"/")[[1]][2]
gmmat_load <- tryCatch(suppressWarnings(library(GMMAT)), error = function(dummy){
lib_folder <- paste0("/",server, "/mbatchou/GMMAT_compute78")
library(GMMAT, lib.loc = lib_folder)
return(0)}
)

# Parameters for simulation
n_total <- 500; n_G <- 11000
F = .01 # fixation index
maf_W1 <- .1
maf_W1_pop1 <- rbeta(1, maf_W1 * (1 - F)/F, (1 - maf_W1) * (1 - F)/F)
maf_W1_pop2 <- rbeta(1, maf_W1 * (1 - F)/F, (1 - maf_W1) * (1 - F)/F)
maf_W2 <- .3
maf_W2_pop1 <- rbeta(1, maf_W2 * (1 - F)/F, (1 - maf_W2) * (1 - F)/F)
maf_W2_pop2 <- rbeta(1, maf_W2 * (1 - F)/F, (1 - maf_W2) * (1 - F)/F)

check_non_poly<- function(G, p_i1, p_i2, n_pop){
	n_total <- n_pop * 2
	tol <- 1e-5
	if (length(G) == 0) return(1)
	if (is.matrix(G)){
		indexes <- (1:dim(G)[2])[(non_poly <- apply(G, 2, sd) == 0 & findInterval(p_i1, c(tol, 1- tol)) == 1 & findInterval(p_i2, c(tol, 1- tol)) == 1 )]
		if (length(indexes) > 0){
			G[1:n_pop,indexes] <- matrix(rbinom(n_pop * length(indexes), 2, rep(p_i1[indexes],n_pop)), n_pop, length(indexes), byrow = TRUE)
			G[(n_pop+1):n_total,indexes] <- matrix(rbinom(n_pop * length(indexes), 2, rep(p_i2[indexes],n_pop)), n_pop, length(indexes), byrow = TRUE)
			G[,indexes] <- check_non_poly(G[,indexes], p_i1[indexes], p_i2[indexes], n_pop)
		} else return(G)
	} else if (sd(G) == 0){
		G <- rbinom(n_total, 2, rep(c(p_i1,p_i2), each = n_pop))
		G <- check_non_poly(G, p_i1, p_i2, n_pop)
	} else return(G)
	return(G)
}

get_G <- function(n_G, n_tot){ # n_G is number of markers (eg 10000)
	p_i <- 	runif(n_G)
	p_i1 <- rbeta(n_G, p_i * (1 - F)/F, (1 - p_i) * (1 - F)/F)
	p_i2 <- rbeta(n_G, p_i * (1 - F)/F, (1 - p_i) * (1 - F)/F)
	G_pop1 <- t(replicate(n_total/2, rbinom(n_G, 2, p_i1)))
	G_pop2 <- t(replicate(n_total/2 + 1, rbinom(n_G, 2, p_i2)))  # Add dummy person for GRM
	G <- rbind(G_pop1, G_pop2)

	# G <- check_non_poly(G, p_i1, p_i2, n_tot / 2)
	# Choose only polymorphic common variants (MAF > .05)
	G <- G[,(apply(G, 2, var) > 0) & (colMeans(G)/2 > 0.05)] 
	paste("There are", n_G - dim(G)[2], "non-polymorphic/rare variants")
	if(dim(G)[2] > 10000) G <- G[,1:10000]
	n_G <<- dim(G)[2]
	return(G)
}

# Add an additional individual and remove him from GRM afterwards
# Compute p_hat as usual but times (maf * n + 1) / (n + 1) --- adding a heterozygote to the count
non_psd <- TRUE
while(non_psd){
	p <- colMeans( (G_mat <- get_G(n_G, n_total)) ) / 2 
	G_mat_c <- scale(G_mat, center = TRUE, scale = sqrt(2*p*(1 - p)))
	GRM <- (tcrossprod(G_mat_c) / n_G )[-(n_total+1),]
	GRM <- GRM[,-(n_total+1)]
	G_mat <- G_mat[-(n_total+1),]
	G_mat_c <- G_mat_c[-(n_total+1),]
	Prec <- tryCatch(chol2inv(chol(GRM)), error = function(dummy) return(FALSE))
	if(is.matrix(Prec)) non_psd <- FALSE
}
Prec <- as(Prec, "sparseMatrix"); rownames(Prec) <- 1:n_total
# tail(eigen(GRM)$va)


# Get pairwise correlations for CARAT
kin <- melt(GRM)[,c(2,1,3)]; kin = kin[kin[,1]<=kin[,2],]

get_sim_pop <- function(n_total, maf_W1, maf_W2, ncase, ncont, int, LG, betas){
	n_sim <- n_total * 10
	mean_W <- (1 - (1-maf_W1)^2) * (1- (1 - maf_W2)^2)
	int <- int - LG * mean_W  # Subtract LG * E(epistatis)
	sim_data = data.frame(sex = rbinom(n_sim, 1, .5),
	age = runif(n_sim, 20, 65),
	Z = rnorm(n_sim),
	W1 = rbinom(n_sim, 2, maf_W1),
	W2 = rbinom(n_sim, 2, maf_W2))
	sim_data$pis <- with(sim_data, plogis(int + LG * ((W1 > 0) & (W2 > 0)) + sex * betas[1] + age * betas[2] + Z * betas[3]))
	sim_data$Y <- rbinom(n_sim, 1, prob = sim_data$pis)
	cases <- sim_data$Y == 1
	sim_data <- rbind(sim_data[cases,][1:ncase,], sim_data[!cases,][1:ncont,])
	return(sim_data)
}

 prop_case1 = .5; var_xb = 11.6;e_lamb = 4.7; prev = .008
 
runme = function (dummy, prop_case1, var_xb, e_lamb, prev)
{
	kt<<- kt + 1; write(kt, "iteration_track.txt")
	ncase1 <- prop_case1 * n_total/2
	ncont1 <- n_total/2 - ncase1
	ncase2 <- n_total/2 - ncase1
	ncont2 <- n_total/2 - ncont1
	
	betas <- sqrt( (var_xb / 3) / c(.25, 168.75, 1)) # (sex, age, Z)
	LG <- log(e_lamb) # increase in penetrance
	interc <- qlogis(prev) - sum(betas * c(.5,42.5,0)) # logit(Preval) - b * E(X) 
	pvals_t1err = NA; pvals_power = NA
	
	# while(any(is.na(pvals_t1err)) | any(is.na(pvals_power))){
	# Sample from population 1 -- no ascertainment done
	sample_data_pop1 <- get_sim_pop(n_total, maf_W1_pop1, maf_W2_pop1, ncase1, ncont1, interc, LG, betas)
	while(any(is.na(sample_data_pop1$Y))){
		sample_data_pop1 <- get_sim_pop(n_total, maf_W1_pop1, maf_W2_pop1, ncase1, ncont1, interc, LG, betas)
	} 
	
	# Sample from population 2 -- no ascertainment done
	sample_data_pop2 <- get_sim_pop(n_total, maf_W1_pop2, maf_W2_pop2, ncase2, ncont2, interc, LG, betas)
	while(any(is.na(sample_data_pop2$Y))){
		sample_data_pop2 <- get_sim_pop(n_total, maf_W1_pop2, maf_W2_pop2, ncase2, ncont2, interc, LG, betas)
	} 
	sample_data <- rbind(sample_data_pop1, sample_data_pop2)
	rm(sample_data_pop1, sample_data_pop2)
	sample_data$idU = 1:nrow(sample_data)	

### MCMC glmm
	prior1 <- list(R = list(V = 1, fix=1),
	G = list(G1 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=50^2)))
	system.time(mMCMC <- MCMCglmm(Y ~ 1 + sex + age + Z, random =~idU, family = "categorical",
	ginverse = list(idU = Prec), prior=prior1,
	data = sample_data, verbose = FALSE, nitt=15e3, thin=100, burnin=5e3))
	pi_hat <- as.vector(predict(mMCMC, marginal = ~idU, type = "response")) # E(pi)

## Generating unassociated Gs
	cols <- (dummy %% 10)*1000 + 1
	Gtot <-  cbind(sample_data$W1, G_mat[,cols:(cols+1000-1)]) # Causal variant + non-causal variants
	n_Gtot <- ncol(Gtot)
	
# Numerator
	cenYres <- sample_data$Y - pi_hat - mean(sample_data$Y - pi_hat)
	T_num_G <- as.vector(crossprod(Gtot, cenYres))^2
	
# Retrospective denominator
	get_denum_retro <- function(G, Matdenum){ 
		mafs <- colMeans(G) / 2
		2 * mafs * (1 - mafs) * Matdenum
	}
	YPhiY <- crossprod(cenYres, GRM) %*% cenYres #marginal mean
	T_denum_G <- get_denum_retro(Gtot, YPhiY)
	pvalstot <- pchisq(T_num_G / T_denum_G, df = 1, lower.tail = FALSE)
	pvals_t1err <- pvalstot[-1]
	pvals_power <- pvalstot[1]

	numend <- sample(1e6,1)
	filend <- paste0(numend, '.txt')
# GMMAT
	mGMMAT <- glmmkin(Y ~ age + sex + Z,  data = sample_data, kins = GRM, family = binomial(link="logit"))
	write.table(cbind(paste0("SNP", 1:(n_Gtot+1)), t(Gtot)), paste0('genofile', filend), row.names=FALSE, col.names=FALSE, quote=F, sep='\t')
	capture.output(glmm.score(mGMMAT, infile = paste0('genofile', filend), outfile = paste0('gmmat_out', filend)), file = "/dev/null")
	if(file.exists(paste0('gmmat_out', filend))){
		pvalstot <- as.numeric(read.table(paste0('gmmat_out', filend), header=TRUE)[,6])
		pvals_t1err <- cbind(pvals_t1err, pvalstot[-1]) 
		pvals_power <- append(pvals_power, pvalstot[1]) 
	} else{
		pvals_t1err <- cbind(pvals_t1err, NA) 
		pvals_power <- append(pvals_power, NA) 
		write(paste('Error in GMMAT', kt), "errorfile.txt", append = TRUE)}

# CARAT
	write.table(cbind(1,1:n_total,2,3,4,sample_data[,c('Y','age','sex','Z')]), paste0('phenofile', filend), sep = '\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
	write(c('Chr','SNP','cm','bp',1:n_total), paste0('genofile', filend), ncolumns = n_total + 4, sep='\t')
	write.table(cbind(1, paste0("SNP", 1:(n_Gtot+1)), 1,50, t(apply(Gtot,2, function(g) unlist(list(c(1,1),c(1,2),c(2,2))[g+1])))),paste0('genofile', filend), row.names=FALSE, col.names=FALSE,  sep='\t', quote = FALSE, append = TRUE)
	write.table(cbind(1,kin),paste0('grmfile', filend),row.names=FALSE, col.names=FALSE,sep=' ')
	com_carat <- paste('./CARAT -p',paste0('phenofile', filend),'-g',paste0('genofile', filend),'-R ',paste0('grmfile', filend),'-o',numend)
	ret <- system(com_carat, ignore.stdout = TRUE, ignore.stderr = TRUE)
	if ((!file.exists(paste0(numend,'_out.txt'))) | (ret !=0)){
		pvals_t1err <- cbind(pvals_t1err, NA) 
		pvals_power <- append(pvals_power, NA) 
		allfile <- list.files()
		file.remove(allfile[grep(numend,allfile)])
		write(paste('Error in CARAT', kt), "errorfile.txt", append = TRUE)
	} else {
		if( length((pvalstot <- as.numeric(read.table(paste0(numend,'_out.txt'), header=TRUE)[,4]))) == 0){
			pvals_t1err <- cbind(pvals_t1err, NA) 
			pvals_power <- append(pvals_power, NA) 
		} else {
			pvals_t1err <- cbind(pvals_t1err, pvalstot[-1]) 
			pvals_power <- append(pvals_power, pvalstot[1]) 
		}
	}
	
# Erase files
		allfile <- list.files()
		file.remove(allfile[grep(numend,allfile)])
	}
	pvalstot <- list(Type1Err = pvals_t1err, Power = pvals_power)
    return(pvalstot)
}

M = seq(.5,.9,.1)
n_M = length(M)
# Will be a list (over Tu/Te) of a list (over rep) of matrices (GMMAT and MCMC)
res = rep(list(list()), n_M) 
n_rep = 100; kt <- 0 # To keep track of iterations

for(i in 1:n_M) {
  res[[i]] = mclapply(1:n_rep, runme, mc.cores = 12, prop_case1 = M[i], var_xb = 11.6, e_lamb = 4.7, prev = .008)
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
col_vec <- rainbow(ncol(power_mat))
pch_vec <- 14 + seq(ncol(power_mat))

pdf("power_plot.pdf", width = 7.5, height = 7)
errbar(props, (y<-power_mat[na_i,1]), y+sds, y-sds, type = "b", pch=pch_vec[1], col=col_vec[1], errbar.col = 'red', ylim = c(0,1), cap = .03,
 ylab = "Power", xlab = "Proportion of variance on logit scale due to polygenic effects/covariates", xaxt = "n" )
for(i in 2:ncol(power_mat)){
	errbar(props, (y<-power_mat[na_i,i]), y+sds, y-sds, type = "b",add=T, cap = .03,pch = pch_vec[i], col=col_vec[i], errbar.col = col_vec[i])
}
axis(1, at = props, label = c("0/100", "20/80","40/60","60/40","80/20","100/0")[na_i])
legend('topright', legend = c('MCMCglmm', 'GMMAT','CARAT'), pch = pch_vec, col = col_vec, cex=1.4, bty="n")
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
method_used = c('MCMCglmm', 'GMMAT','CARAT')
power_mat = matrix(NA, 6, 3)
type1err <- matrix(NA, 6, 3, dimnames = list(NULL, method_used))

index <- (1:6)[lapply(res,length)>0]
for(k in index){
assign(paste0('p_vals',k,'_t1err'), vector("list", 3)) # List for 3 methods
assign(paste0('p_vals',k,'_pow'), vector("list", 3)) # List for 3 methods
}

for(i in 1:3){
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
n_pvals <- NULL
for(k in 1:length(index)){
	eval(parse(text=paste0('type1err[',index[k],',] <- unlist(lapply(p_vals',index[k],'_t1err, pvals_analysis, al = .05))')))
	eval(parse(text=paste0('n_pvals[',k,'] <- length(p_vals',index[k],'_t1err[[1]])')))
}
type1err; t(sapply(n_pvals, function(x).05 + c(-1,1)*2*sqrt(.05*(1-.05)/x)))
png("t1err_plot%01d.png", width = 12, height = 9, units = 'in', res=200)
par(mfrow = c(2,2))
for(k in 1:length(index)){
	eval(parse(text=paste0('lapply(1:3, draw_pval, p_vals',index[k],'_t1err, method_used)')))
	title(paste('Polygenic effect / Covariates:',props[index[k]],'/',100-props[index[k]]), outer=TRUE, cex=1.1, line = -1.5)
}
dev.off()


################################## for 6 settings of Tu/Tot
n_p=rep(0, n_M)
method_used = c('MCMCglmm', 'GMMAT','CARAT')
power_mat <- matrix(NA, n_M, 3)
type1err <- matrix(NA, n_M, 3, dimnames = list(NULL, method_used))

for(k in 1:n_M){
assign(paste0('p_vals',k,'_t1err'), vector("list", 3)) # List for 3 methods
assign(paste0('p_vals',k,'_pow'), vector("list", 3)) # List for 3 methods
}

for(i in 1:3){
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
n_pvals <- NULL
for(k in 1:n_M){
	eval(parse(text=paste0('type1err[',k,',] <- unlist(lapply(p_vals',k,'_t1err, pvals_analysis, al = .05))')))
	eval(parse(text=paste0('n_pvals[',k,'] <- length(p_vals',k,'_t1err[[1]])')))
}
write.table(type1err, row.names=FALSE, quote=FALSE); t(sapply(n_pvals, function(x).05 + c(-1,1)*2*sqrt(.05*(1-.05)/x)))
png("t1err_plot%d.png", width = 12, height = 9, units = 'in', res=200)
par(mfrow = c(2,2))
for(k in 1:n_M){
	eval(parse(text=paste0('lapply(1:3, draw_pval, p_vals',k,'_t1err, method_used)')))
	title(paste('Polygenic effect / Covariates:',props[k],'/',100-props[k]), outer=TRUE, cex=1.1, line = -1.5)
}
dev.off()
