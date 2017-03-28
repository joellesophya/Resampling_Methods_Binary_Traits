split.direct.sum<-function(x){  

	if(is.na(x)){
		return(NULL)
	}else{
		x <- gsub("\n", "", x, fixed=TRUE)  # required if the random formula is very long and broken over lines
		openB<-gregexpr("\\(", x)[[1]]
		closeB<-gregexpr("\\)", x)[[1]]
		true_openB<-openB
		true_closeB<-closeB

		for(i in 1:length(openB)){
			dist<-outer(openB, closeB, function(x,y){y-x})
			dist[which(dist<0)]<-Inf
			dist<-which(dist==min(dist), arr.ind=T)[1,] 
			true_openB[i]<-openB[dist[1]]
			true_closeB[i]<-closeB[dist[2]]
			openB<-openB[-dist[1]]
			closeB<-closeB[-dist[2]]
		}

		plus<-gregexpr("\\+", x)[[1]]
		internals<-matrix(mapply(function(x,y){(plus>x & plus<y)}, x=true_openB, y=true_closeB), length(plus), length(true_openB))
		rterms<-strsplit(x, "")[[1]]
		rterms[plus[which(rowSums(internals)!=0)]]<-"leaveMCMCleave"
		rterms<-paste(rterms, collapse="")
		rterms<-strsplit(rterms, " *\\+ *")[[1]]
		rterms<-gsub("leaveMCMCleave", "+", rterms)
		return(rterms)
	}
}

predict_MCMC <- function (object, newdata = NULL, marginal = object$Random$formula, 
type = "response", interval = "none", level = 0.95, it = NULL, 
posterior = "all", verbose = FALSE, ...) 
{
	rm.obs <- c()
	rcomponents <- split.direct.sum(as.character(object$Random$formula)[2])
	mcomponents <- split.direct.sum(as.character(marginal)[2])
	marginalise <- rep(as.numeric(rcomponents %in% mcomponents), 
	object$Random$nrt)
	if (posterior == "mean") {
		object$VCV <- matrix(colMeans(object$VCV), 1, ncol(object$VCV))
		object$Sol <- matrix(colMeans(object$Sol), 1, ncol(object$Sol))
		it <- 1
	}
	if (posterior == "mode") {
		object$VCV <- matrix(posterior.mode(object$VCV, ...), 1, ncol(object$VCV))
		object$Sol <- matrix(posterior.mode(object$Sol, ...), 1, ncol(object$Sol))
		it <- 1
	}
	if (is.null(it)) {
		if (posterior == "all") {
			it <- 1:nrow(object$Sol)
		}
	}
	object$Sol <- object$Sol[it, , drop = FALSE]
	object$VCV <- object$VCV[it, , drop = FALSE]
	if (is.null(object$Random$nfl) == FALSE) {
		st <- c(1, cumsum(object$Random$nrl * object$Random$nfl) + 
		1)
		st <- st[-length(st)]
		end <- cumsum(object$Random$nrl * object$Random$nfl)
		keep <- unlist(mapply(st[which(marginalise == 0)], 
		end[which(marginalise == 0)], FUN = ":"))
	}
	else {
		keep <- NULL
	}
	object$Sol <- object$Sol[, c(1:object$Fixed$nfl, object$Fixed$nfl + 
	keep), drop = FALSE]
	W <- cBind(object$X, object$Z)
	W <- W[, c(1:object$Fixed$nfl, object$Fixed$nfl + keep), 
	drop = FALSE]
	post.pred <- t(apply(object$Sol, 1, function(x) {
		(W %*% x)@x    # Gets the x slot from vector (W%*%x) which is from dgCMatrix class so post pred is regular matrix
	}))
	if (type == "response") {
		if (any(object$family != "gaussian" & object$family != 
					"cengaussian")) {
			post.var <- buildV(object, marginal = marginal, 
			diag = TRUE, it = NULL, posterior = "all", 
			...)
		}
		super.trait <- 1:length(object$Residual$family)
		cnt <- 1
		i <- 1
		while (i <= length(super.trait)) {
			if (grepl("multinomial", object$Residual$family[i])) {
				nm <- as.numeric(substr(object$Residual$family[i], 
				12, nchar(object$Residual$family[i]))) - 
				1
				super.trait[i + 1:nm - 1] <- cnt
				i <- i + nm
				cnt <- cnt + 1
			}
		}
		normal.logistic <- function(mu, v) {
			int.foo <- function(x, mu, v) {
				plogis(x) * dnorm(x, mu, sqrt(v))
			}
			integrate(int.foo, qnorm(0.0001, mu, sqrt(v)), 
			qnorm(0.9999, mu, sqrt(v)), mu, v)[[1]]
		}
		normal.multilogistic <- function(mu, v) {
			int.foo <- function(x, mu, v, i) {
				(exp(x[i])/(1 + sum(exp(x)))) * prod(dnorm(x, 
				mu, sqrt(v)))
			}
			res <- 1:length(mu)
			for (i in 1:length(mu)) {
				res[i] <- cubature::adaptIntegrate(int.foo, 
				qnorm(0.0001, mu, sqrt(v)), qnorm(0.9999, 
				mu, sqrt(v)), mu = mu, v = v, i = i)[[1]]
			}
			res
		}
		for (k in unique(super.trait)) {
			if (any(grepl("multinomial", object$Residual$family))) {
				keep <- which(object$error.term %in% which(super.trait == 
				k))
				size <- as.numeric(substr(object$family[keep], 
				12, nchar(object$family[keep])))
				for (j in 1:nrow(post.pred)) {
					prob <- matrix(post.pred[j, keep], length(keep)/sum(super.trait == 
					k), sum(super.trait == k))
					pvar <- matrix(post.var[j, keep], length(keep)/sum(super.trait == 
					k), sum(super.trait == k))
					post.pred[j, keep] <- t(sapply(1:nrow(prob), 
					function(x) {
						normal.multilogistic(prob[x, ], pvar[x, 
						])
					})) * size
				}
			}
		}
	}
	if (length(rm.obs) > 0) {
		post.pred <- post.pred[, -rm.obs]
	}
	if (is.matrix(post.pred)) {
		pred <- matrix(colMeans(post.pred), dim(post.pred)[2], 
		1)
	}
	else {
		pred <- matrix(post.pred, length(post.pred), 1)
	}
	rownames(pred) <- 1:dim(pred)[1]
	return(pred)
}
