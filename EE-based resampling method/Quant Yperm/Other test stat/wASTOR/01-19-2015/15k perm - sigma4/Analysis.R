setwd("C:/Users/JOELLE/Google Drive/Reasearch Mary Sara _ SAMSUNG//C Code//Project/Permutation based methods/Binary trait/Results/")

ASTOR_t=read.table('01-19-2015/EE method - log - 15k perm - sigma4/ASTOR_statistics.txt',header=F)
ASTOR_t=unlist(c(ASTOR_t));names(ASTOR_t)=c()
p=read.table('01-19-2015/EE method - log - 15k perm - sigma4/sim_pvals.txt',header=F)
p=unlist(c(p));names(p)=c()

######################
## Analysis
summary(ASTOR_t)
p=1-pchisq(ASTOR_t,1)

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

## Plot for one-tailed test
uni2= rank(p_vals,ties.method='max'); names(uni2)=c()
plot( -log10(uni2/n2),-log10(p_vals), pch=1, ylab='-log(Observed p-value)', xlab='-log(Expected p-value)',
      cex=.7)
abline(a=0,b=1,col="red")
title("Scenario 5: EE method with logistic",line=3)
title(expression(paste(sigma[a]^2,"=4 and 15000 permutation replicates")),line=2)

##########################
## Combining both on plot
######################
plot( -log10(uni/n),-log10(p), pch=1, ylab='-log(Observed p-value)', xlab='-log(Expected p-value)',
      main='Logistic model using LMM method with ASTOR
      (asymptotic and perm based p-values) ', cex=.7,  col='blue')
points(-log10(uni2/n2),-log10(p_vals), pch=1 ,cex=.7,col='black')
abline(a=0,b=1,col="red")

