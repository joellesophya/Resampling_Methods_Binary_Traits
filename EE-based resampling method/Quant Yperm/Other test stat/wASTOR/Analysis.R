
ASTOR_t=read.table('ASTOR_statistics.txt',header=F)
ASTOR_t=unlist(c(ASTOR_t));names(ASTOR_t)=c()
p=read.table('sim_pvals.txt',header=F)
p=unlist(c(p));names(p)=c()

######################
## Analysis
summary(ASTOR_t)
pASTOR=1-pchisq(ASTOR_t,1)

get_type1err <- function(p){
	# summary(p)
	n=length(p)
	
	hist(p, main="Histogram of permutation-based p-values",freq = F,breaks=100)
	abline(h=1,col="red")

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

get_type1err(p); get_type1err(pASTOR)
png("pvals_plot%01d.png", width = 6, height = 4.5, units = 'in', res=150)
plot_p(p, "EE method"); plot_p(pASTOR, "ASTOR")
dev.off()
