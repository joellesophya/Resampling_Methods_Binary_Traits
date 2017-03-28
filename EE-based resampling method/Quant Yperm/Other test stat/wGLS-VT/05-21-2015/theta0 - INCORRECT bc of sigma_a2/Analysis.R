# EE, LOG, 15k-25k Y perm, sigma_a^2=0, using VT fam
setwd("C:/Users/JOELLE/Google Drive/Reasearch Mary Sara _ SAMSUNG//C Code//Project/Permutation based methods/Binary trait")

p_vals=read.table('Results/05-21-2015/EE - VT fam - log - theta0/pvals_25000.txt',header=F)
p_vals=unlist(c(p_vals));names(p_vals)=c()

########################### For p_vals from permutation method
summary(p_vals)
n=length(p_vals);n
hist(p_vals, main="Histogram of permutation-based p-values",freq = F,breaks=100)
title(expression(paste(theta[u]," = 0")),line=.7)
abline(h=1,col="red")

alpha=c(.005,.01,.05)
result=data.frame("alpha"=alpha, "Err_rate"=c(0,0,0),"SE"=c(0,0,0),"p_value"=c(0,0,0))
for(i in 1:3){
  a=alpha[i]
  result$p_value[i]=prop.test(x = sum(p_vals<a),n,p = a)$'p.val'
  result$Err_rate[i]=round(as.numeric(prop.test(x = sum(p_vals<a),n,p = a)$'est'),4)
  result$SE[i]=round(sqrt(mean(p_vals<a)*(1-mean(p_vals<a))/n),4)
};result


## Which should be the minimum p_hat to get .05 in CI
s=function(p){p+qt(1-.05/2,n-1)*sqrt(p*(1-p)/n)-.05}
uniroot(s,c(0,1))$root

## Plot for one-tailed test
uni2= rank(p_vals,ties.method='max'); names(uni2)=c()
plot( -log10(uni2/(n+1)),-log10(p_vals), type="n", ylab='Observed (-log10P)', xlab='Expected (-log10P)')
title("Plot of observed vs expected -log10 p-values",line=3)
title(expression(paste(theta[u],"=0 and 25k permutation replicates")),line=2)

a=seq(0.00001,1,by=.0001)
low=-log10(a+qnorm(1-.05/2)*sqrt(a*(1-a)/n))
high=-log10(a-qnorm(1-.05/2)*sqrt(a*(1-a)/n))
polygon(c(-log10(a), rev(-log10(a))), c(high, rev(low)), col ='gray', border = NA)
points(-log10(uni2/(n+1)),-log10(p_vals), pch=16,cex=.3)
abline(a=0,b=1,col="red")
