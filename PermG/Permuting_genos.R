library(mvtnorm)
library(parallel)
library(reshape)

# Parameters for simulation
N_fam <- 45
n_person <- 22
n_total <- N_fam * n_person
n_marker <- 30
n_perm <- 1e3
MAF <- runif(n_marker, .1,.3) # population allele frequencies
tol <- 1e-6

# Make pedigree
ped1 <- data.frame(fam = NA, ind = 1:n_total, fa = NA, mo = NA)
fa_ind <- c(NA,NA,NA,1,1,NA,1,NA,1,NA,4,4,4,6,6,6,8,8,10,10,10,10)
mo_ind <- c(NA,NA,NA,2,2,NA,2,NA,2,NA,3,3,3,5,5,5,7,7,9,9,9,9)
ped1[,c('mo','fa')] <- cbind(rep(mo_ind, N_fam)+ n_person*rep(0:(N_fam-1), each = n_person),
rep(fa_ind, N_fam)+ n_person*rep(0:(N_fam-1), each = n_person))
ped1[is.na(ped1)] <- 0
ped1$fam <- rep(1:N_fam, each = n_person)
# write.table(ped1,'pedigree',row.names=FALSE, col.names=FALSE,sep='\t')

# Generate data
phi <- as.matrix(read.table("add_phi_mat_n22")); colnames(phi) <- 1:n_person
# Phimat <- diag(N_fam) %x% phi
# kin <- melt(diag(N_fam)%x%phi)[,c(2,1,3)]; kin = kin[kin[,1]<=kin[,2],]
# Get fam for row and col indID
# kin$famrow <- ((kin[,1]-1)%/%n_person) + 1; kin$famcol <- ((kin[,2]-1)%/%n_person) + 1
# kin <- kin[kin$famrow==kin$famcol,]
# phi_inv <- solve(phi) #22x22 matrix
sex <- rep(c(1,2,2,1,2,1,2,1,2,1,2,1,1,2,1,2,1,2,2,1,2,1),N_fam) #constant

# Get chol, V
chol_phi <- chol(phi)
C_inv <- solve(chol_phi)
W <-  colSums(C_inv)
WtW_inv <- 1/(N_fam * sum(W^2))
V <- svd(diag(n_total) - tcrossprod(rep(W,N_fam)) * WtW_inv)
keep_V <- V$d>tol  # sum(keep_V) == (n_person - 1)
V1 <- V$u[,keep_V]

peddrop <- function (parentmat, n_G, MAF) { # parentmat is matrix of parental info
	nfound <- sum(is_found <- ((parentmat[, 1] == 0) & (parentmat[, 2] == 0)))
	rec_f <- runif(n_G - 1, 0.05, 0.15)
	npeople <- dim(parentmat)[1]
	haplos <- matrix(0, npeople, 2)
	haplos[is_found,] <- matrix(1:(2 * nfound), nfound, 2, byrow=T)
	haplos <- rep(list(haplos), n_G)
	
	for (i in 1:npeople) {
		#/* When both parents are in the pedigree drop alleles down the pedigree */
		if ((parentmat[i, 1] != 0) & (parentmat[i, 2] != 0)) {
			#/* Randomly choose chrom. 0 or 1 to transmit to kid at each locus
			trchrom <- sample(1:2, 2, replace = T)
			#// For first GT haplotype
			haplos[[1]][i,] <- sapply(1:2, function(a) haplos[[1]][parentmat[i, a], trchrom[a]])
			# For next GT haplotype
			for(i_G in 2:n_G){
				# Transmit from same chrom. 1-rf frac of the time 
				k <- runif(2) < rec_f[i_G - 1]
				trchrom <- sapply(1:2, function(a){
					if(k[a] == TRUE) setdiff(1:2, trchrom[a])
					else trchrom[a]
				})
				haplos[[i_G]][i,] <- sapply(1:2, function(a) haplos[[i_G]][parentmat[i, a], trchrom[a]])
			} 
		}
	}
	state <- lapply(1:n_G, function(a) as.numeric(runif(2 * nfound) < MAF[a])) #1 if minor, 0 if major allele
	genos <- sapply(1:n_G, function(a) rowSums(matrix(state[[a]][haplos[[a]][,]], npeople, 2)))/2 # MAC/2 for each person at each marker (n_person x n_marker)
	return(genos)
}

compute_statistic <- function(G, Y){
	## Choose a statistic to use 
	lm(Y~G)$coef[2]
}
runme = function (dummy)
{
	kt<<- kt + 1; write(kt, "iteration_track.txt")
	# Generate Y under some model
	a_var <- 4
	U <- c(t(rmvnorm(N_fam, sigma = a_var * phi)))
	Y <- 3 + 2 * sex + U
	
	## Generate the Gs and the permuted genotypes
	Gmat<- do.call(rbind, lapply(1:N_fam, function(a) peddrop( ped1[1:n_person, c("fa","mo")], n_marker, MAF ))) # N x n_marker matrix of tightly linked markers [rf U(.05-.15)]
	Gmat_stat <- apply(Gmat, 2, compute_statistic, Y=Y)
	p_hat <- sapply(1:n_marker, function(a) sum(crossprod(W, matrix(Gmat[,a], n_person, N_fam)))*WtW_inv)
	Z <- sapply(1:n_marker, function(a) crossprod(V1, crossprod(diag(N_fam)%x% C_inv, scale(Gmat, center = p_hat, scale = FALSE)[,a])))
	
	Gperm_stat <- t(sapply(1:n_perm, function(a){
		perm <- sample(1:nrow(Z), replace = F)
		G_perm <- scale(crossprod(diag(N_fam)%x% chol_phi, V1 %*% Z[perm,]), center = -p_hat, scale = FALSE) # N x n_marker matrix
		# Expensive to store the whole Gperm for all n_perm so perhaps compute statistic directly in here
		apply(G_perm, 2, compute_statistic, Y=Y)
	}))
	est_pval <- (rowSums(t(Gperm_stat) > Gmat_stat) + 1) / (n_perm + 1)


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
