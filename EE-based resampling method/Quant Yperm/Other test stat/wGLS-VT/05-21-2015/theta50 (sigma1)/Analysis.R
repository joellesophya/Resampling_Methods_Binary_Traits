# 100 Y rep, 15k-35k Y perm, sigma_a^2=4, using VT statistic with RVs
setwd("C:/Users/JOELLE/Google Drive/Reasearch Mary Sara _ SAMSUNG//C Code//Project/Permutation based methods/Binary trait")

p_vals=read.table('Results/05-21-2015/EE - VT fam - log - theta50/pvals_25000.txt',header=F)
p_vals=unlist(c(p_vals));names(p_vals)=c()

########################### For p_vals from permutation method
summary(p_vals)
n2=length(p_vals)
hist(p_vals, main="Histogram of permutation-based p-values",freq = F,breaks=100)
abline(h=1,col="red")

alpha=c(.005,.01,.05)
result=data.frame("alpha"=alpha, "Err_rate"=c(0,0,0),"p_value"=c(0,0,0),"SE"=c(0,0,0))
for(i in 1:3){
  a=alpha[i]
  result$p_value[i]=prop.test(x = sum(p_vals<a),n2,p = a)$'p.val'
  result$Err_rate[i]=round(as.numeric(prop.test(x = sum(p_vals<a),n2,p = a)$'est'),4)
  
  result$SE[i]=round(sqrt(mean(p_vals<a)*(1-mean(p_vals<a))/n2),4)
};result

## Which should be the minimum p_hat to get .05 in CI
s=function(p){p+qt(1-.05/2,n-1)*sqrt(p*(1-p)/n)-.05}
uniroot(s,c(0,1))$root

## Plot for one-tailed test
uni2= rank(p_vals,ties.method='max'); names(uni2)=c()
plot( -log10(uni2/n2),-log10(p_vals), pch=1, ylab='-log(Observed p-value)', xlab='-log(Expected p-value)',
      cex=.7)
abline(a=0,b=1,col="red")
title("Plot of observed vs expected -log p-values",line=3)
title(expression(paste(sigma[a]^2,"=4 and 15k permutation replicates")),line=2)

##########################
## Combining both on plot
######################
plot( -log10(uni/n),-log10(p), pch=1, ylab='-log(Observed p-value)', xlab='-log(Expected p-value)',
      main='Logistic model using LMM method with ASTOR
      (asymptotic and perm based p-values) ', cex=.7,  col='blue')
points(-log10(uni2/n2),-log10(p_vals), pch=1 ,cex=.7,col='black')
abline(a=0,b=1,col="red")
