######################
## Analysis

get_type1err <- function(p){
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
	result
}

## QQ plot of the p-values
plot_p <- function(p_vals, method){
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

n_pvals <- 25 * 1e3
main_folder <- "Scenario 4"
i <- 1; hnum <- ifelse(i==0, '', "h10")
# main_folder <- paste0(main_folder, hnum)
pvals_folder <- letters[seq(2)]

pvals <- lapply(pvals_folder, function(let){
	if(length(pvals_folder) ==1){
		p <- read.table(paste0(main_folder,"/","pvals_", n_pvals,".txt"),header=F)
		p=unlist(c(p));names(p)=c()
		png(paste0(main_folder,"/pvals_plot_",n_pvals,".png"), width = 6, height = 4.5, units = 'in', res=150)
		plot_p(p, paste(n_pvals, "p-values"))
		dev.off()
		get_type1err(p)
	} else {
	# subdir <- paste0(main_folder,"/", main_folder, let, hnum,"/")
	subdir <- paste0(main_folder,"/", main_folder, hnum, let,"/")
	p <- read.table(paste0(subdir,"pvals_", n_pvals,".txt"),header=F)
	return(p)
	}
})

if(length(pvals_folder) !=1){
	ptot <- unlist(pvals)
	png(paste0(main_folder,"/pvals_plot", hnum,"_",n_pvals,".png"), width = 6, height = 4.5, units = 'in', res=150)
	plot_p(ptot, paste(n_pvals, "p-values"))
	dev.off()
	get_type1err(ptot)
}



##### Dist of naff in Ypi vs in Y
#ytot <- cbind(Y, mu_y + tcrossprod(pre_mult_mat, t(replicate(n_perm, sample(delt)))))
 ytot <- cbind(Y, (replicate(n_perm, sample(Y)))) #naive perm
 y <- apply(ytot, 2, function(y) tapply(y, ped1$fam, sum))
 yavg <- rowMeans(y[,-1])
 
 hist(yavg,col=rgb(0.8,0.8,0.8,0.5), xlim=c(-.5,.5), ylim=c(0,7), freq=F)
 hist(y[,1], col=rgb(0.1,0.1,0.1,0.5), freq=F, add = T)

 png("Dist_naff_in_Ypi_naive_perm.png")
 par(mfrow=c(2,3))
 sapply(floor(seq(1,45,length = 6)), function(nf) {
	hist(y[nf,-1], freq = F, main = paste("Family", nf), xlim = c(-1, 16),
	xlab ="Number of affected in Y replicate")
	legend(x=9.4,y=.13, legend="Obs. Y", lwd=1, col='red', bty='n', cex=.85, seg.len = 1)
	abline(v = y[nf,1], col = "red")
	})
 dev.off()
 
