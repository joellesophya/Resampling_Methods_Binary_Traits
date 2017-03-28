library(mvtnorm)
library(parallel)
library(reshape)

library(MCMCglmm)

library(mailR)
sendmemail <- function(){
to <- "joelle.mbatchou@gmail.com"
subject <- c("Batch job done!")
body <- tail(strsplit(getwd(),"/")[[1]],1); body <- paste("Folder: ", body)
mailControl=list(host.name = "aspmx.l.google.com", port = 25)
send.mail(from=to,to=to,subject=subject, body = body, smtp=mailControl,
authenticate = F, send = T)
}

# Parameters for simulation
N_fam <- 45
n_person <- 22
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
kin <- melt(diag(N_fam)%x%phi)[,c(2,1,3)]; kin = kin[kin[,1]<=kin[,2],]
# Get fam for row and col indID
kin$famrow <- ((kin[,1]-1)%/%n_person) + 1; kin$famcol <- ((kin[,2]-1)%/%n_person) + 1
kin <- kin[kin$famrow==kin$famcol,]
write.table(kin[,c(4,1:3)],'kinship',row.names=FALSE, col.names=FALSE,sep=' ')
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
  S_XB <- sqrt(Tot * (1 - Tu) / (363.67+.25+1)) # Ratio of target variance to var(X)
  beta <- c(.05, 1, .5) * S_XB
  LG <- log(1.1) # increase of 10% in the odds of being a case
  int <- qlogis(0.1) - beta[1] * 31.86 -  beta[2] * 1.5 - LG * .1425  # Prevalence - b * E(X) - lambda * .1425
  
  age <- rep(c(73,75,46,43,40,46,40,43,47,51,18,21,15,15,12,9,13,17,24,21,18,14),N_fam) + runif(n_total,-1.5,1.5)
  Z <- rnorm(n_total)
  U <- c(t(rmvnorm(N_fam, sigma = add_var * phi))) 
  G1 <- rowSums(peddrop( ped1[,c("fa","mo")], 0.1 ))
  G2 <- rowSums(peddrop( ped1[,c("fa","mo")], 0.5 ))
  
  Y <- rbinom(n_total, size = 1, prob = plogis(int + U + LG * ((G1 > 0) & (G2 > 0)) +	age*beta[1] + sex*beta[2] + Z*beta[3]))
  fam_data = data.frame(Y=Y, age = age, sex = sex, Z = Z, idU = 1:n_total)
  
  ### MCMC glmm
  # MCMCtime <- proc.time()
  prior1 <- list(R = list(V = 1, fix=1),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=10^2)))
  mMCMC<-MCMCglmm(Y ~ 1 + age + sex + Z, random =~idU, family = "categorical",
                  ginverse = list(idU = Prec), prior=prior1,
                  data = fam_data, verbose = F, nitt=dummy * 250e3, thin=500, burnin=10e3)
  
  k <- ((16*sqrt(3))/(15*pi))^2
  est_add_var <- posterior.mode(mMCMC$VCV[,1] / (1 + k))
  est_beta <- posterior.mode(mMCMC$Sol[,1:4] / sqrt(1 + k))
  
  ### Can't do check of (Y - pi) - Prec * u* = 0 b/c u* is not computed in MCMCglmm at each iteration 
  # --> increase number of MCMC iter so that the posterior for u is well estimated ???
  # con_pi <- as.vector(predict(mMCMC, posterior = "mode", type = "response")) # conditional on posterior mode of parameters (beta, u) -- E(pi|u)
  mar1 <- as.vector(predict(mMCMC, marginal = ~idU, type = "response")) # E(pi)
  ## Get marginal pi by approx from formula 2.14 in mcmcglmm course notes
  # ru <- rowSums(mMCMC$VCV)
  # X <- as.matrix(mMCMC$X)
  # b <- mMCMC$Sol
  # xb <- tcrossprod(b, X)
  # mar2<- colMeans(plogis(xb - .5 * ru * tanh(xb *  (1 + 2 * exp(-0.5 * ru))/6)))   #  RMSE = sqrt(sum((mar1-mar2)^2))  -- closer
  
  # Generate 1000 G
  m <- 5000
  G_sub <- replicate(m, rowSums(peddrop( ped1[,c("fa","mo")], 0.3 ))) # n * 1000 matrix
  
  # Numerator
  T_num <- as.vector(crossprod(G_sub, Y - mar1))^2
  # T_num1 <- as.vector(crossprod(G_sub, Y - con_pi))^2
  # for(i in 1:2) assign(paste("T_num",i+1,sep=""), eval(parse(text = paste("as.vector(crossprod(G_sub, Y - mar",i,"))^2", sep=""))))
  
  # Prospective Denum
  # get_denum1 <- function(pivec){  # Based on taylor approx of score and mean pi -- works for all
    # Gamm <- diag(pivec * (1 - pivec))
    # Matcen <- Gamm - crossprod(Gamm, X) %*% solve(crossprod(X, crossprod(Gamm, X))) %*% crossprod(X, Gamm)
    # colSums(G_sub * (Matcen %*% G_sub))
  # }
  # T_denum1 <- get_denum1(con_pi) 
  # for(i in 1:2) assign(paste("T_denum",i+1,sep=""), eval(parse(text = paste("get_denum1(mar",i,")", sep=""))))
  
  # Create fisher info matrix then invert and get element corresponding to GT effect (n_cov+1)
  # get_denum2 <- function(pivec){  # Based on including dependency of u* on snp effect 
  # Gamm <- diag(pivec * (1 - pivec))
  # sig_hat <- posterior.mode(mMCMC$VCV[,1])
  # Sig_inv <- diag(N_fam) %x% phi_inv / sig_hat
  # W_mat <- Sig_inv + Gamm
  # C_mat <- crossprod(2 * diag(n_total) - W_mat^2, crossprod(W_mat, Gamm))
  # Hbb <- crossprod(X,Gamm) %*% (C_mat - diag(n_total)) %*% X
  # Hbg <- crossprod(X,Gamm) %*% (C_mat - diag(n_total)) %*% G_sub
  # Hbs <- - crossprod(X, t(W_mat)) %*% crossprod(Sig_inv, U)
  # Hgg <- colSums(G_sub * (crossprod(Gamm, C_mat - diag(n_total)) %*% G_sub))
  # Hgs <- - c(crossprod(G_sub, t(W_mat)) %*% crossprod(Sig_inv, U))
  # Hss <- 1/sig_hat^2 * crossprod(U,Sig_inv) %*% (crossprod(4 * diag(n_total) - W_mat^2, crossprod(W_mat, Sig_inv)) - diag(n_total)) %*% U -
  # 2/sig_hat^2 * crossprod(U,Sig_inv) %*% crossprod(W_mat, Y - pivec) +
  # 1/sig_hat^2 * crossprod(U,Sig_inv) %*% (crossprod(W_mat, Sig_inv) %*% W_mat - Sig_inv) %*% (Y - pivec - crossprod(Sig_inv, U))
  
  # I_mat <- lapply(split(seq(m),seq(m)), function(x){
  # -rbind(cbind(Hbb, Hbg[,x], Hbs),
  # c(Hbg[,x], Hgg[x], Hgs[x]),
  # c(Hbs, Hgs[x], Hss))
  # })
  
  
  # }
  # T_denumV1 <- get_denum2(con_pi)
  # for(i in 1:2) assign(paste("T_denumV",i+1,sep=""), eval(parse(text = paste("get_denum2(mar",i,")", sep=""))))
  
  # Retrospective version
  Phimat <- diag(N_fam) %x% phi
  get_denum_retro <- function(G, Matdenum){ 
    maf <- mean(G) / 2 #(sum(G) + 1) / ((n_total + 1) * 2)  --> Not necessary since G is common variant
    2 * maf * (1 - maf) * Matdenum
  }
  YPhiY <- crossprod(Y - mar1, Phimat) %*% (Y - mar1) #marginal mean
  Tr_denumR <- apply(G_sub, 2, get_denum_retro, YPhiY)
  
  #### QQ plots for comparison
  # draw_quant <- function(T){
  # y <- qchisq(ppoints(length(T_num1)), df = 1)
  # qqplot(y, T, xlab = bquote("Expected" ~ chi[1]^2 ~ "quantiles"), ylab = "Observed quantiles")
  # abline(a = 0, b = 1, col = "blue") #;qqline(T, distribution = function(p) qchisq(p, df = 1), col = 2)
  # }
  # # Prospective
  # i=1
  # eval(parse(text = paste("draw_quant(T_num",i,"/T_denum",i,")",sep=""))) # based on taylor approx
  # eval(parse(text = paste("draw_quant(T_num",i,"/T_denumV",i,")",sep=""))) # based on fisher info
  # Retrospective
  # draw_quant(T_num2 / Tr_denumR)
  
  #### P-values
  draw_pval <- function(p_vals){
    n = length(p_vals)
    uni2= rank(p_vals, ties.method='max'); names(uni2)=c()
    plot(-log10(uni2/(n+1)),-log10(p_vals), type="n", ylab=expression('Observed (-log'[10]*' p-value)'), xlab=expression('Expected (-log'[10]*' p-value)'))
    a=1:n
	high <- qbeta(0.025, a, rev(a))
	low <- qbeta(0.975, a, rev(a))
	polygon(-log10(c(a/n,rev(a/n))), -log10(c(high, rev(low))), col ='gray', border = NA)
	points(-log10(uni2/(n+1)),-log10(p_vals), pch=16,cex=.3) 
    title(expression(bold('Plot of observed vs. expected -log'[10]*' p-values')),line=3.2)
    title(paste("Setting:", Tu,"/", Tot, "- #MCMC iter:", 250*dummy,"k"), line=1.5)
    abline(a=0,b=1,col="red")
  }
  pvals_analysis <- function(p_vals){
    n <- length(p_vals); i <- 1
    alpha <- c(.005,.01,.05)
    result <- data.frame("alpha" = alpha, "Err_rate" = NA,"SE" = NA,"p_value" = NA)
    for(a in alpha){
      result$p_value[i] <- prop.test(x = sum(p_vals < a), n, p = a)$'p.val'
      result$Err_rate[i] <- round(as.numeric(prop.test(x = sum(p_vals < a), n, p = a)$'est'),4)
      result$SE[i] <- round(sqrt(mean(p_vals < a)*(1-mean(p_vals < a)) / n),4)
      i <- i + 1
    }; result
  }
  
  # Prospective
  # for(i in 1:3) assign(paste("pvals",i,sep=""), eval(parse(text = paste("pchisq(T_num",i,"/ T_denum",i,", df = 1, lower.tail = F)", sep=""))))
  # for(i in 1:3) assign(paste("pvalsV",i,sep=""), eval(parse(text = paste("pchisq(T_num",i,"/ T_denumV",i,", df = 1, lower.tail = F)", sep=""))))
  # Retrospective
  pvalsR<- pchisq(T_num / Tr_denumR, df = 1, lower.tail = F)
  
  pdf(paste("pvals_Retro",dummy,"_",Tu*100,".pdf"))
  # for(i in 1:3){
    # eval(parse(text = paste("draw_pval(pvals",i,",",i,")", sep="")))
    # eval(parse(text = paste("pvals_analysis(pvals",i,")", sep="")))
  # }
  draw_pval(pvalsR)
  result_out = pvals_analysis(pvalsR)
  dev.off()
  
  return(result_out)
}

########################### 
################## Analysis
M = rbind(c(0,1.1),c(.2,1.2),c(.4,1.2), c(.6,1.3),c(.8,1.3),c(1,1.3))
n_M = dim(M)[1]
# Will be a list (over Tu/Te) of a list (over rep) of matrices (GMMAT and MCMC)
res = rep(list(list()), n_M) 
n_rep = 2; kt <- 0 # To keep track of iterations

for(i in 1:n_M) {
  res[[i]] = mclapply(1:n_rep, runme, mc.cores = 2, Tu= M[i,1], Tot = M[i,2])
	save.image() #Save the res object
}

# sendmemail()

# p_valsMCMC <- vector("list", n_M) # List of size n_M
# est_varMCMC <- vector("list", n_M)
# m <- length(res[[1]][[1]]) # Number of pvals + est add var

# for(i in 1:n_M) {
# for(j in 1:n_rep){
# p_valsMCMC[[i]] <- append(p_valsMCMC[[i]], res[[i]][[j]][-m])
# est_varMCMC[[i]] <- append(est_varMCMC[[i]], res[[i]][[j]][m])
# }
# }

# pdf("Summaries of pvals_10pct.pdf")

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
# main = "Variance estimates over 200 Y realizations",
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

