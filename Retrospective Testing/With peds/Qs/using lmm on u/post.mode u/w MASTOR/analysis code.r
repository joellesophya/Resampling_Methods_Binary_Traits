########################### 
################## Analysis
plot_power <- function(power_mat, index, meth_list){
sds <- 2*sqrt(.25 / n_p)
n_method <- length(meth_list)
props <- ((0:10)*.1)[index]

pdf("power_plot.pdf")
errbar(props, (y<-power_mat[,1]), y+sds, y-sds, type = "b", pch=16, col='red', errbar.col = 'red', ylim = c(0,1), cap = .03,
 ylab = "Power", xlab = "Proportion of variance on logit scale due to polygenic effects/covariates", xaxt = "n" )
for(i in 2:n_method){
	errbar(props, (y<-power_mat[,i]), y+sds, y-sds, type = "b",add=T, cap = .03,pch = c(15,17,18)[i-1], col=c('blue','lightgreen','purple')[i-1], errbar.col = c('blue','lightgreen','purple')[i-1])
}
axis(1, at = props, label = paste(props*100, (1-props)*100, sep = "/") ) 
legend('topright', horiz = F, legend = meth_list, pch = c(16,15,17,18)[1:n_method] , col = c('red','blue','lightgreen','purple')[1:n_method], cex=1.4, bty="n")
dev.off()
}

pvals_analysis <- function(p_vals, al){
    n <- length(p_vals); i <- 1
    alpha <- c(.005,.01,.05)
    result <- data.frame("alpha" = alpha, "Err_rate" = NA,"SE" = NA,"p_value" = NA)
    for(a in alpha){
      result$p_value[i] <- prop.test(x = sum(p_vals < a), n, p = a)$'p.val'
      result$Err_rate[i] <- round(as.numeric(prop.test(x = sum(p_vals < a), n, p = a)$'est'),4)
      result$SE[i] <- round(sqrt(mean(p_vals < a)*(1-mean(p_vals < a)) / n),4)
      i <- i + 1
    }; result$Err_rate[alpha==al]
  }

draw_pval <- function(ind, p_vals_list, meth_list, percent){
  p_vals <- p_vals_list[[ind]]
  method <- meth_list[ind]
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

########################### 
########################### 
n_ind <- length(index <- ((1:11)[round((0:10)*.1, 1) %in% M[,1]])[sapply(res,length) > 0])

n_p=rep(0, n_ind)
n_meth <- length( (method_used = c('MCMCglmm_BlueGLS','CERAMIC')) )
power_mat = matrix(NA, n_ind, n_meth+1)
type1err <- matrix(NA, n_ind, n_meth, dimnames = list(NULL, method_used))
n_pvals <- NULL

a_var_est <- matrix(NA, n_rep, length(index))
for(k in index){
assign(paste0('p_vals',k,'_t1err'), vector("list", n_meth)) # List for methods
assign(paste0('p_vals',k,'_pow'), vector("list", n_meth + 1)) # List for methods
}


for(i in 1:(n_meth + 1)){
	for(j in 1:n_rep){
		for(k in 1:length(index)){
			if(is.double(unlist(res[[k]][[j]]))){
				if(i <= n_meth) eval(parse(text=paste0('p_vals',index[k],'_t1err[[',i,']] <- append(p_vals',index[k],'_t1err[[',i,']] , res[[',k,']][[',j,']]$Ty[,',i,'])')))
				eval(parse(text=paste0('p_vals',index[k],'_pow[[',i,']] <- append(p_vals',index[k],'_pow[[',i,']] , res[[',k,']][[',j,']]$Po[',i,'])')))
				if (i ==1){
					n_p[k] = n_p[k] + 1
					a_var_est[j,k] <- res[[k]][[j]]$add
				}
			}
		}
	}
}

# index <- 1:6
level <- .05
# Power
for(k in 1:length(index)){
	eval(parse(text=paste0('power_mat[',k,',] <- unlist(lapply(p_vals',index[k],'_pow, function(x) mean(x<level)))')))
}
library(Hmisc)
plot_power(power_mat, index, method_used)

# Type 1 Error
props <- seq(0,100,by=10)
for(k in 1:length(index)){
	eval(parse(text=paste0('type1err[',k,',] <- unlist(lapply(p_vals',index[k],'_t1err, pvals_analysis, al = .05))')))
	eval(parse(text=paste0('n_pvals[',k,'] <- length(p_vals',index[k],'_t1err[[1]])')))
}
type1err;t(sapply(n_pvals, function(x).05 + c(-1,1)*2*sqrt(.05*(1-.05)/x)))
png("t1err_plot%01d.png", width = 12, height = 9, units = 'in', res=200)
par(mfrow = c(2,1))
for(k in 1:length(index)){
	eval(parse(text=paste0('lapply(1:n_meth, draw_pval, p_vals',index[k],'_t1err, method_used)')))
	title(paste('Polygenic effect / Covariates:',props[index[k]],'/',100-props[index[k]]), outer=TRUE, cex=1.1, line = -1.5)
}
dev.off()
