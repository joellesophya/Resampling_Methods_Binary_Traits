# 100 Y rep
# 10000 Y perm
# sigma_a^2=1 using all 500 markers with EE-logistic model

ASTOR_t=read.table('ASTOR_statistics.txt',header=F)
ASTOR_t=unlist(c(ASTOR_t));names(ASTOR_t)=c()
p=read.table('sim_pvals.txt',header=F)
p=unlist(c(p));names(p)=c()

######################
## Analysis
summary(ASTOR_t)
pASTOR=1-pchisq(ASTOR_t,1)

get_type1err <- function(p){
	summary(p)
	n=length(p)
	hist(p, main="Histogram of permutation-based p-values",freq = F,breaks=100)
	abline(h=1,col="red")

	alpha=c(.005,.01,.05)
	result=data.frame("alpha"=alpha, "Err_rate"=c(0,0,0),"p_value"=c(0,0,0),"SE"=c(0,0,0))
	for(i in 1:3){
		a=alpha[i]
		result$p_value[i]=prop.test(x = sum(p<a),n,p = a)$'p.val'
		result$Err_rate[i]=round(as.numeric(prop.test(x = sum(p<a),n,p = a)$'est'),4)

		result$SE[i]=round(sqrt(mean(p<a)*(1-mean(p<a))/n),4)
	};result
}
> get_type1err(p); get_type1err(pASTOR)
  alpha Err_rate   p_value     SE
1 0.005   0.0052 0.5055743 0.0003
2 0.010   0.0108 0.0687073 0.0005
3 0.050   0.0509 0.3720708 0.0010
  alpha Err_rate    p_value     SE
1 0.005   0.0052 0.46591051 0.0003
2 0.010   0.0108 0.06214129 0.0005
3 0.050   0.0511 0.28135641 0.0010

## Plot for one-tailed test
plot_p <- function(p){
	uni= rank(p,ties.method='max'); names(uni)=c()
	plot( -log10(uni/length(p)),-log10(p), pch=1, ylab='-log(Observed p value)', xlab='-log(Expected p-value)',
	main='Asymptotic p-values from logistic model using ASTOR
	(empirical type I error rate of ) ', cex=.7)
	abline(a=0,b=1,col="red")
}

## Which should be the minimum p_hat to get .05 in CI
s=function(p){p+qt(1-.05/2,n-1)*sqrt(p*(1-p)/n)-.05}
uniroot(s,c(0,1))$root

##########################
## Combining both on plot
######################
plot( -log10(uni/n),-log10(p), pch=1, ylab='-log(Observed p-value)', xlab='-log(Expected p-value)',
      main='Logistic model using LMM method with ASTOR
      (asymptotic and perm based p-values) ', cex=.7,  col='blue')
points(-log10(uni2/n2),-log10(p_vals), pch=1 ,cex=.7,col='black')
abline(a=0,b=1,col="red")
