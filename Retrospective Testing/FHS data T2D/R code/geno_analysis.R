# install.packages("data.table", repos = "https://Rdatatable.github.io/data.table", type = "source")
pckg_list <- c("MCMCglmm", "parallel", "data.table", "Matrix") 
dummy <- lapply(pckg_list, require, character.only = TRUE)
load("MCMC_null_res.RData")

ctphi <- t(chol(as.matrix(phi)))
setnames(kin <- data.table(melt(as.matrix(phi))), 1:3, c("person2","person1","kincoef"))
kin <- kin[person1<=person2 & (pheno$unipedno[person1] == pheno$unipedno[person2])]
kin[, fam := pheno$unipedno[person1]]
# For CERAMIC -- same for all realizations
kin[kincoef == 1, kincoef:= 0]; kin[, kincoef := kincoef/2]
write.table(kin[,.(fam, person1, person2, kincoef)],'grmfile', row.names=F, col.names=F,sep=' ')

system("mkdir pvalsout")

##### Get genotypes
shareid_retained = pheno$shareid
shareid_geno = fread("genodata/ShareID_allchr.txt")$V1   # File about column ID for genotype files
ind_retained = match(shareid_retained,shareid_geno)    # Get column location in genotype file for each individual
# Get file that has SNP information as well as file that has QC results for all SNPs
gene_info <- cbind(fread("SNPs_name_mapping.txt", col.names = c("SNP_name","SNP_rs","CHR","BP"), colClasses = "character", select = 1:4), fread("xwretainsnps_map.txt",col.names = "keep"))
# Remove SNPs that do not pass QC criterions for call rate, Mendelian error rate and MAF.
gene_info <- gene_info[keep== 1]; gene_info[,keep:=NULL]
# Remove SMPs where p-value testing association with cohort ID <=1e-7
xwWQLS  <- fread("WQLSpval.txt", col.names=c("SNP","pval"))
gene_info <- gene_info[gene_info$SNP_name %in% xwWQLS[ pval>1e-7, SNP]]
rm(xwWQLS)
# Remove SNPs in X chromosome or with no chromosome information
gene_info <- gene_info[CHR != "X" & CHR !="---"]
# Convert column classes
gene_info[, `:=`(CHR = as.integer(CHR), BP = as.numeric(BP))]

# Function to analyze chromosome i
chr_analysis <- function(i){
	system(paste0("mkdir CERAMIC_chr", i))  # To accomodate runnning CERAMIC in parallel
# Files to use for raw genotype data, CERAMIC input genos and p-values 
	input_file_name = paste0("./genodata/genodata_chr",i,".txt")  # obtained from ./Framingham Files/work_sheng
	CERAMIC_filename = paste0("./CERAMIC_chr",i,"/genofile")
	pvals_out = paste0("./pvalsout/genoP_chr",i,".txt")

# Read genotype file
	gendata <- fread(input_file_name, select = ind_retained)
	gendata <- gendata[, paste0("V",ind_retained), with=F] # re-order columns according to order in data set
	name_snps <- fread(paste0("./genodata/SNPID_chr",i,".txt"), header = F, select = 1)$V1  # File with SNP labels
# Filter based on QC results
	gendata <- gendata[name_snps %in% gene_info[, SNP_name]]
	name_snps <- name_snps[name_snps %in% gene_info[, SNP_name]]
	gendata[gendata == -4] <- -9 # convert missing values for CERAMIC
	fwrite(gendata, file = CERAMIC_filename, col.names = F, row.names = F, sep='\t') # For CERAMIC input genofile
	fwrite(list(name_snps) , file = paste0("./CERAMIC_chr",i,"/snpfile"), col.names = F, row.names = F, sep='\t') # CERAMIC SNP ID file

# Prepare genofile for retrospective test
	gendata[gendata == -9] <- NA_real_ # convert missing values
	gendata <- as(t(gendata), "sparseMatrix")  # columns are the SNPs
# Copy gene matrix and replace missing values by 0
	gendata_non0 <- gendata
	gendata_non0[is.na(gendata)] <- 0

####### RETRO_MCMC 
	t0 <- proc.time()  ## How long it takes to analyze each chromosome
# create phenotypic residual for each marker
	big_pheno <- as(replicate(ncol(gendata), M2),"sparseMatrix")
	big_pheno[is.na(gendata)] <- 0  # Remove phenotypic info for individuals with missing GT for each marker
	big_pheno <- scale(big_pheno, center = T, scale = F)  # Mean center the residual -- E(G^t(Y-mu))=0 under retro model

# test numerator
	N2 <- colSums(gendata_non0 * big_pheno)^2
# test denomirator
	GLS_sig <- diag(crossprod(gendata_non0, P) %*% gendata_non0)/(nrow(P) - 1)  # GLS estimate of variance for Gs
	D2 <- GLS_sig * rowSums(crossprod(big_pheno, ctphi)^2) # YtPhiY = sum(CY)^2 for C = chol(Phi)
	pvalsM <- pchisq(N2 / D2, df = 1, lower.tail = F)
	t1 <- proc.time()
	
# Run CERAMIC from terminal (each chromosome has its own folder)
	write.table( pheno[,.(unipedno, .I, pheno, sex, BMI_ave)], paste0('./CERAMIC_chr',i,'/','phenofile'), sep = '\t', quote=F, col.names=F, row.names=F)
	com_ceramic <- paste0('(cd ./CERAMIC_chr',i,' && exec ../BATMAN -p phenofile -g genofile -k ../grmfile -n snpfile)')
	ret <- system(com_ceramic, ignore.stdout = T, ignore.stderr = T)
	pvalsC <- fread(paste0('./CERAMIC_chr',i,'/BATMANtest.pvalues'), skip = 1, col.names = c("b","SNP","CERAMIC"))[,-1]
	pvalsC[, RETRO_MCMC := pvalsM]
	fwrite(pvalsC, pvals_out, sep = "\t")
	if(i==22) system(paste0('mv ./CERAMIC_chr',i,'/BATMANtest.phenoestimates ./')) # Keep results from fitting null model (same for all chr)
	system(paste0('rm -R CERAMIC_chr', i))
	return(t1-t0) 
}

res = mclapply(1:22, chr_analysis, mc.cores = 7)
save.image()

######## Analyze results
pvals <- do.call("rbind", lapply(1:22, function(i){
	pvals_out = paste0("./pvalsout/genoP_chr",i,".txt")
	pmat <- fread(pvals_out)
	pmat[,CHR := i]  # Add chromosome information
}))
nrow(pvals)
ind_in_geno_info <- match(pvals$SNP, gene_info$SNP_name)  # Recover SNP info (name and BP position) from gene_info 
pvals[, rs := gene_info[ind_in_geno_info, SNP_rs]]
pvals[, BP := gene_info[ind_in_geno_info, BP]]
setkey(pvals, RETRO_MCMC) # order by ascending p-values

# Get top 10 and 20 p-values (could make a function to automatically get top k p-values)
top20 <- pvals[1:20,.(rs, CHR,  BP,RETRO_MCMC, CERAMIC)]
setnames(top20, c("SNP", "Chr","Position","RETRO_MCMC", "CERAMIC"))
setkey(top20, Chr)
top20[,minp := c("RETRO","CER")[as.numeric(RETRO_MCMC>=CERAMIC)+1]] # Which of the two has smallest p-value
top20[, `:=`(RETRO_MCMC = format(RETRO_MCMC, scientific = T, digits=2), CERAMIC = format(CERAMIC, scientific = T, digits=2))]
x<-sapply(1:20, function(i) humarray::nearest.gene(top20[i, Chr], top20[i, Position],1,side="either")) # get nearest gene info
top20[, Nearest_Gene := x]
top20 <- top20[,.(SNP, Chr,Position,Nearest_Gene,RETRO_MCMC, CERAMIC, minp)]
fwrite(top20, "top20.txt", sep="\t")

top10 <- pvals[1:10,.(rs, CHR,  BP,RETRO_MCMC, CERAMIC)]
setnames(top10, c("SNP", "Chr","Position","RETRO_MCMC", "CERAMIC"))
setkey(top10, Chr)
top10[,minp := c("RETRO","CER")[as.numeric(RETRO_MCMC>=CERAMIC)+1]]
top10[, `:=`(RETRO_MCMC = format(RETRO_MCMC, scientific = T, digits=2), CERAMIC = format(CERAMIC, scientific = T, digits=2))]
x<-sapply(1:10, function(i) humarray::nearest.gene(top10[i, Chr], top10[i, Position],1,side="either"))
top10[, Nearest_Gene := x]
top10 <- top10[,.(SNP, Chr,Position,Nearest_Gene,RETRO_MCMC, CERAMIC, minp)]
fwrite(top10, "top10.txt", sep="\t", quote=F)

####### Some plots
## QQ Plots with confidence bands (my own function)
plot_QQ(pvals$RETRO_MCMC, "RETRO-MCMC")
plot_QQ(pvals$CERAMIC, "CERAMIC")

library(qqman)
manhattan(pvals, chr = "CHR", bp = "BP", p = "RETRO_MCMC", snp = "rs", main = "With RETRO-MCMC")
manhattan(pvals, chr = "CHR", bp = "BP", p = "CERAMIC", snp = "rs", main = "With CERAMIC")

## Genomic inflation factor
chisq <- qchisq(1-pvals$RETRO_MCMC,1) #.99
chisq <- qchisq(1-pvals$CERAMIC,1) #1.00
(lamb <- median(chisq)/qchisq(0.5,1))
1 + (lamb - 1) * (1/577 + 1/1190)/(1/1000+1/1000)
