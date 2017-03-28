# Binary Y no cov &asc with WSS
n_rv = 150
setwd("C:/Users/JOELLE/Google Drive/Reasearch Mary Sara _ SAMSUNG//C Code//Project/Permutation based methods/Binary trait")

p_vals1=read.table('Results/09-30-2015/Scenario 6 h60/pvals_20000.txt',header=F)[-(1:2),]
p_vals2=read.table('Results/10-01-2015/Scenario 6 h60a/pvals_20000.txt',header=F)
p_vals3=read.table('Results/10-01-2015/Scenario 6 h60b/pvals_20000.txt',header=F)
p_vals4=read.table('Results/10-04-2015/Scenario 6 h60a/pvals_20000.txt',header=F)
p_vals5=read.table('Results/10-04-2015/Scenario 6 h60b/pvals_20000.txt',header=F)
p_vals6=read.table('Results/10-04-2015/Scenario 6 h60c/pvals_20000.txt',header=F)
p_vals1=unlist(c(p_vals1));names(p_vals1)=c()
p_vals2=unlist(c(p_vals2));names(p_vals2)=c()
p_vals3=unlist(c(p_vals3));names(p_vals3)=c()
p_vals4=unlist(c(p_vals4));names(p_vals4)=c()
p_vals5=unlist(c(p_vals5));names(p_vals5)=c()
p_vals6=unlist(c(p_vals6));names(p_vals6)=c()
p_vals=c(p_vals1,p_vals2,p_vals3,p_vals4,p_vals5,p_vals6); rm(p_vals1,p_vals2,p_vals3,p_vals4,p_vals5,p_vals6)

########################### For p_vals from permutation method
summary(p_vals)
n=length(p_vals);n/n_rv
hist(p_vals, main="Scenario 6 - with Weighted Sum Stat",freq = F,breaks=100, xlab="p value")
title(expression(paste(" Proportion of variance = 60%")),line=.7)
abline(h=1,col="red")

alpha=c(.005,.01,.05)
result=data.frame("alpha"=alpha, "Err_rate"=c(0,0,0),"SE"=c(0,0,0),"p_value"=c(0,0,0))
for(i in 1:3){
  a=alpha[i]
  result$p_value[i]=prop.test(x = sum(p_vals<a),n,p = a)$'p.val'
  result$Err_rate[i]=round(as.numeric(prop.test(x = sum(p_vals<a),n,p = a)$'est'),4)
  result$SE[i]=round(sqrt(mean(p_vals<a)*(1-mean(p_vals<a))/n),4)
};result


## Which should be the minimum p_hat to get al in CI
al=.05
s_l=function(p){p+qt(1-al/2,n-1)*sqrt(p*(1-p)/n)-al}
s_u=function(p){p-qt(1-al/2,n-1)*sqrt(p*(1-p)/n)-al}
c(uniroot(s_l,c(0,1))$root, uniroot(s_u,c(0,1))$root)

## Plot for one-tailed test
uni2= rank(p_vals,ties.method='max'); names(uni2)=c()
plot( -log10(uni2/(n+1)),-log10(p_vals), type="n", ylab=expression('Observed (-log'[10]*' p-value)'), xlab=expression('Expected (-log'[10]*' p-value)'))
title(expression(bold('Plot of observed vs. expected -log'[10]*' p-values')),line=3.2)
title(expression(paste("Scenario 4 with WST")),line=2)

a=seq(0.00001,1,by=.0001)
low=-log10(a+qnorm(1-.05/2)*sqrt(a*(1-a)/n))
high=-log10(a-qnorm(1-.05/2)*sqrt(a*(1-a)/n))
polygon(c(-log10(a), rev(-log10(a))), c(high, rev(low)), col ='gray', border = NA)
points(-log10(uni2/(n+1)),-log10(p_vals), pch=16,cex=.3)
abline(a=0,b=1,col="red")
