# setwd("~/Framingham Files/work_Joelle")
pckg_list <- c("MCMCglmm", "parallel", "data.table", "Matrix") 
# "matrixStats",, ,"reshape")
dummy <- lapply(pckg_list, require, character.only = TRUE)

## Phenotype and covariates file
pheno <- fread("phenocov.txt")
# Remove people with missing phenotype or covariates
## 625 cases with non-missing covariates, and 1413 controls with non-missing covariates = 2038 total
pheno <- pheno[pheno != 0 & !is.na(BMI_ave)]

### MCMC chain has bad mixing when data contains many pedigrees with all people being all cases or all controls
# --> For such pedigrees, keep only one individual
# --> If using marginal mean as estimate of pi, since it only depends on add_var_est, could keep just one individual per such pedigree and fit the null then estimate pi for all individuals using X, beta_hat and add_var_hat
ped_rem <- pheno[, sd(pheno),pedno]
ped_rem <- ped_rem[V1 == 0, pedno] # pedigrees with all cases/cont
# pheno[pedno%in%ped_rem, .N, pedno] # Count of people in the peds (2-8 inds)
# Create variable for ID w/in pedigree
pheno[, IDped:= 1:.N, by = pedno]
pheno <- pheno[!(pedno %in% ped_rem & IDped > 1)]

# Create unique pedno variabl eand reorder pheno by ascending pedno
n_person_fam <- pheno[,.N,pedno]
setorder(n_person_fam, pedno)
unaff <- n_person_fam[pedno==0,N]
ord <- order(pheno$pedno)
pheno <- rbind(pheno[ord[-(1:unaff)]],pheno[ord[1:unaff]]) 
pheno[, unipedno := c(rep(n_person_fam[pedno!=0,.I], n_person_fam[pedno!=0,N]), n_person_fam[pedno!=0,.N] + 1:unaff)]

# Kinship matrix
phimat <- fread("kinfile.txt", col.names = c("pedno","ind1","ind2","kin"))
phimat <- phimat[ind1 %in% pheno$shareid & ind2 %in% pheno$shareid]
phimat[, sind1 := match(ind1, pheno$shareid)]
phimat[, sind2 := match(ind2, pheno$shareid)]
# Get inverse of relatedness matrix (i.e. 2 * kinship coeff)
phi <- as(as.matrix(dcast(phimat[,], sind1~sind2, value.var = "kin", fill = 0)[,-1, with = F]), "sparseMatrix")
phi <-  2 * (phi + t(phi) )
diag(phi) <- 1
Prec <- solve(phi)
P <- Prec - outer(rowSums(Prec), colSums(Prec)) / sum(Prec) # For retro test
rownames(Prec) <- unique(phimat$sind1); colnames(Prec) <- unique(phimat$sind2)
rm(phimat)

# Get variables of interest
pheno_sub <- pheno[,.(pheno, sex, BMI_ave)]
pheno_sub[, `:=`(pheno = pheno - 1, unishareid = .I)]


# Fit null model
prior1 <- list(R = list(V = 1, fix=1),
G = list(G1 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=10^2)))
tMCMC <- system.time(mMCMC <- MCMCglmm(pheno ~ 1 + sex + BMI_ave, random =~unishareid, family = "categorical", ginverse = list(unishareid = Prec), prior=prior1, data = pheno_sub, verbose = F, nitt=500e3, thin=1000, burnin=10e3, pr = TRUE))
ncov <- mMCMC$Fixed$nfl

# Estimates from MCMCglmm
ck <- ((16*sqrt(3))/(15*pi))^2
effS <- effectiveSize(mMCMC$VCV)[1]
est_betas <- as.matrix(mMCMC$Sol[,1:ncov]) /sqrt(1 + ck)#niter*ncov
est_u <- as.matrix(mMCMC$Sol[,-(1:ncov)]) /sqrt(1 + ck)
est_s2 <- drop(mMCMC$VC[,1])/(1 + ck)

M2 <- pheno_sub$pheno - plogis(drop(mMCMC$X %*% colMeans(est_betas) + colMeans(est_u)))

rm(Prec)
save.image("FHS_gwas_results.RData")

# Done with prospective model
pheno[,.N,pedno][pedno==0 | N == 1, sum(N)] # 802 unaffected
pheno[,.N,pedno][!(pedno==0 | N == 1),sum(N)] # 965 related
pheno[,.N,pedno][!(pedno==0 | N == 1),.N] # 179 pedigrees
 pheno[,.N,pheno] # 1190 controls and 577 cases

