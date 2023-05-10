
rm(list = ls())

library(data.table)
library(lmtest)
library(parallel)
library(tidyr)
# library(latex2exp)

n.train <- 125
n.test <- 10000
k <- 45
esigma <- 150
reps <- 1000
cores <- 192

y.scale <- sqrt(sum((1:k)^2)+esigma^2)

set.seed(04251978)

gen.true.sigma2 <- function(j){
	sigma2 <- (1/y.scale)^2*(sum(((k-j):0)^2) + esigma^2)
	return(sigma2)
}

true.sigma2 <- sapply(1:k,gen.true.sigma2)

append.df <- function(df.list,dataframe){
	df.list <- append(df.list,list(copy(dataframe)))
	return(df.list)
}

m.copies.df <- function(dataframe,m){
	df.list <- list(copy(dataframe))
	for(i in 2:m){
		df.list <- append.df(df.list,dataframe)
	}
	return(rbindlist(df.list))
}

gen.data <- function(n,k,esigma){
	x <- rnorm(k*n)
	x <- matrix(x,nrow=n)
	x <- data.table(x)
	kmat <- data.table(matrix(rep(rev(1:k),each=n),nrow=n))
	x <- x*kmat
	xnames <- paste0("x.",1:k)
	setnames(x,xnames)
	x[,e:=esigma*rnorm(n)]
	x[,y:=(1/y.scale)*rowSums(.SD),.SDcols=c(xnames,"e")]
	return(x)
}

gen.nonstoch <- function(j,train.data,test.data){
	copies <- floor(nrow(test.data)/nrow(train.data))
	used.x <- m.copies.df(train.data[,1:j],copies)
	alt.data <- cbind(used.x,test.data[,(j+1):ncol(test.data)])
	xnames <- paste0("x.",1:k)
	alt.data[,y:=(1/y.scale)*rowSums(.SD),.SDcols=c(xnames,"e")]
	return(alt.data)
}

formula.j <- function(j){
	paste("y",paste0("x.",1:j,collapse=' + '),sep=' ~ ')
}

train.data <- function(j,dataset){
	eqn <- formula.j(j)
	reg <- lm(eqn,data=dataset)
	return(reg)
}

gen.aux.reg <- function(reg){
	dataset <- data.table(reg$model)
	e <- as.numeric(reg$residuals)
	hat <- as.numeric(lm.influence(reg)$hat)
	sigma2 <- (1+hat)*e^2/(1-hat)
	dataset[,y:=sigma2]
	eqn <- formula.j(ncol(dataset)-1)
	aux <- lm(eqn,data=dataset)
	return(aux)
}

gen.true.resid <- function(j,dataset){
	xvars <- names(dataset)[1:j]
	yhat <- dataset[,(1/y.scale)*rowSums(.SD),.SDcols=xvars]
	e <- dataset[,y]-yhat
	return(e)
}

gen.sample.sigma2 <- function(j,dataset){
	e <- gen.true.resid(j,dataset)
	e2 <- e^2
	sigma2 <- sum(e2)/nrow(dataset)
	return(sigma2)
}

gen.est.sigma2 <- function(reg,dataset){
	yhat <- as.numeric(predict(reg,dataset))
	ybar <- mean(dataset[,y])
	e <- dataset[,y]-yhat
	e2 <- e^2
	sigma2 <- sum(e2)/nrow(dataset)
	return(sigma2)
}

gen.hat.sigma2 <- function(reg,type){
	e <- as.numeric(reg$residuals)
	hat <- as.numeric(lm.influence(reg)$hat)
	if(type=="lin"){
		sigma2 <- sum(e^2/(1-hat))/length(e)
	} else if(type=="square"){
		sigma2 <- ((reg$df.residual)/length(e))*sum(e^2/(1-hat)^2)/length(e)
	} else if(type=="delta"){
		delta <- length(e)*hat/ncol(reg$model)
		delta[delta>4] <- 4
		sigma2 <- sum(e^2/(1-hat)^delta)/length(e)
	}
	return(sigma2)
}

calc.r2 <- function(reg,dataset){
	yhat <- as.numeric(predict(reg,dataset))
	ybar <- mean(dataset[,y])
	e <- dataset[,y]-yhat
	rss <- sum(e^2)
	tss <- sum((dataset[,y]-ybar)^2)
	r2 <- 1-rss/tss
	return(r2)
}

calc.nonstoch.r2 <- function(reg){
	yhat <- as.numeric(reg$fitted.values)
	ybar <- mean(yhat)
	e <- as.numeric(reg$residuals)
	hat <- as.numeric(lm.influence(reg)$hat)
	rss <- sum((1+hat)*e^2/(1-hat))
	tss <- sum((yhat+e-ybar)^2)
	r2 <- 1-rss/tss
	return(r2)
}

calc.stoch.r2 <- function(reg,dataset){

  xmatrix <- data.table(reg$model)
  y <- xmatrix[,y]
  xnames <- copy(names(xmatrix))
  xnames[1] <- "const"
  setnames(xmatrix,xnames)
  xmatrix[,const:=1]
  xmatrix <- as.matrix(xmatrix)
  xpxinv <- solve(t(xmatrix) %*% xmatrix)

  xO <- dataset[,xnames[2:length(xnames)],with=FALSE]
  xO[,const:=1]
  setcolorder(xO,xnames)
  xO <- as.matrix(xO)
  hO <- xO %*% xpxinv
  hO <- hO %*% t(xmatrix)

  hat <- as.numeric(lm.influence(reg)$hat)
  e <- as.numeric(reg$residuals)

  mean.y <- mean(y)
  sq <- as.matrix((y-mean.y)^2)
  sqs <- as.numeric(hO %*% sq)

  sigma2 <- e^2/(1-hat)
  sigma2 <- as.matrix(sigma2)
  sigma2hat <- as.numeric(hO %*% sigma2)
  h2sigma2hat <- as.numeric(hO^2 %*% sigma2)

  ess <- sum(sqs)
  rss <- sum(sigma2hat)
  hss <- sum(h2sigma2hat)

  r2 <- 1-(rss+hss)/ess

  return(r2)
}

predict.loo <- function(i,eqn,dataset){
	obs.list <- 1:nrow(dataset)
	reg <- lm(eqn,data=dataset[obs.list!=i])
	yhat <- predict(reg,dataset[obs.list==i])
	return(yhat)
}

calc.loo.r2 <- function(j,dataset){
	eqn <- formula.j(j)	
	yhat <- as.numeric(sapply(1:nrow(dataset),predict.loo,eqn,dataset))
	ybar <- mean(dataset[,y])
	e <- dataset[,y]-yhat
    rss <- sum(e^2)
    tss <- sum((dataset[,y]-ybar)^2)
    r2 <- 1-rss/tss
    return(r2)
}

calc.mape <- function(reg,dataset,type){

	true.e2 <- true.sigma2[reg$rank-1]
	est.e2 <- as.numeric(reg$residuals)^2

	if(type=="sample"){
		e2 <- gen.true.resid(reg$rank-1,dataset)^2
	} else if(type=="est"){
		e2 <- est.e2
	} else if(type=="adj"){
		e2 <- length(true.e2)*est.e2/reg$df.residual
	} else {
		hat <- as.numeric(lm.influence(reg)$hat)
		if(type=="hat"){
			e2 <- est.e2/(1-hat)
		} else if(type=="hat2"){
			e2 <- (reg$df.residual/length(true.e2))*est.e2/(1-hat)^2
		} else if(type=="delta"){
			delta <- length(true.e2)*hat/reg$rank
			delta[delta>4] <- 4
			e2 <- est.e2/(1-hat)^delta
		}
	}

	mape <- mean(abs(e2-true.e2)/true.e2)

	return(mape)

}

gen.hat.sigma2 <- function(reg,type){
	e <- as.numeric(reg$residuals)
	hat <- as.numeric(lm.influence(reg)$hat)
	if(type=="lin"){
		sigma2 <- sum(e^2/(1-hat))/length(e)
	} else if(type=="square"){
		sigma2 <- ((reg$df.residual)/length(e))*sum(e^2/(1-hat)^2)/length(e)
	} else if(type=="delta"){
		delta <- length(e)*hat/ncol(reg$model)
		delta[delta>4] <- 4
		sigma2 <- sum(e^2/(1-hat)^delta)/length(e)
	}
	return(sigma2)
}

run.simulation <- function(iteration){
	cat(paste0("rep ",iteration,"\n"))
	training.data <- gen.data(n.train,k,esigma)
	test.data <- gen.data(n.test,k,esigma)
	nonstoch.data <- lapply(1:k,gen.nonstoch,training.data,test.data)

	train.sigma2 <- sapply(1:k,gen.sample.sigma2,training.data)
	test.sigma2 <- sapply(1:k,gen.sample.sigma2,test.data)

	k.regs <- lapply(1:k,train.data,training.data)
	train.estsigma2 <- sapply(k.regs,gen.est.sigma2,training.data)
	nonstoch.estsigma2 <- mapply(gen.est.sigma2,k.regs,nonstoch.data)
	test.estsigma2 <- sapply(k.regs,gen.est.sigma2,test.data)

	df.resid <- as.numeric(sapply(k.regs,`[`,"df.residual"))
	train.adjsigma2 <- train.estsigma2*n.train/df.resid

	train.hatsigma2 <- sapply(k.regs,gen.hat.sigma2,"lin")
	train.hat2sigma2 <- sapply(k.regs,gen.hat.sigma2,"square")
	train.deltasigma2 <- sapply(k.regs,gen.hat.sigma2,"delta")

	sigmas <- data.table(i=iteration,k=1:k,true.sigma2=true.sigma2,train.sigma2=train.sigma2,test.sigma2=test.sigma2,train.estsigma2=train.estsigma2,
		nonstoch.estsigma2=nonstoch.estsigma2,test.estsigma2=test.estsigma2,train.adjsigma2=train.adjsigma2,train.hatsigma2=train.hatsigma2,
		train.hat2sigma2=train.hat2sigma2,train.deltasigma2=train.deltasigma2)

	train.r2s <- sapply(k.regs,calc.r2,training.data)
	nonstoch.r2s <- mapply(calc.r2,k.regs,nonstoch.data)
	test.r2s <- sapply(k.regs,calc.r2,test.data)

	nonstoch.hetero.r2s <- sapply(k.regs,calc.nonstoch.r2)
	stoch.nonstoch.r2s <- mapply(calc.stoch.r2,k.regs,nonstoch.data)
	oos.r2s <- sapply(k.regs,calc.stoch.r2,test.data)
	# loo.r2s <- sapply(1:k,calc.loo.r2,training.data)

	r2s <- data.table(i=iteration,k=1:k,train=train.r2s,nonstoch=nonstoch.r2s,test=test.r2s,
		nonstoch.hetero=nonstoch.hetero.r2s,stoch.nonstoch=stoch.nonstoch.r2s,oos=oos.r2s)

	mape.sample <- sapply(k.regs,calc.mape,training.data,"sample")
	mape.est <- sapply(k.regs,calc.mape,training.data,"est")
	mape.adj <- sapply(k.regs,calc.mape,training.data,"adj")
	mape.hat <- sapply(k.regs,calc.mape,training.data,"hat")
	mape.hat2 <- sapply(k.regs,calc.mape,training.data,"hat2")
	mape.delta <- sapply(k.regs,calc.mape,training.data,"delta")

	mapes <- data.table(i=iteration,k=1:k,sample=mape.sample,est=mape.est,adj=mape.adj,hat=mape.hat,hat2=mape.hat2,delta=mape.delta)

	simulation.data <- list("sigmas"=sigmas,"r2s"=r2s, "mapes"=mapes)

	cat(paste0("completed rep ",iteration,"\n"))
	
	return(simulation.data)

}

sims <- mclapply(1:reps,run.simulation,mc.cores=cores)
saveRDS(sims,'sims.homosked45_125.RDS')
sigmas <- lapply(sims,`[[`,"sigmas")
sigmas <- rbindlist(sigmas)

sigmas[,sample.ape:=abs(true.sigma2-train.sigma2)/true.sigma2]
sigmas[,unadj.ape:=abs(true.sigma2-train.estsigma2)/true.sigma2]
sigmas[,adj.ape:=abs(true.sigma2-train.adjsigma2)/true.sigma2]
sigmas[,hat.ape:=abs(true.sigma2-train.hatsigma2)/true.sigma2]
sigmas[,hat2.ape:=abs(true.sigma2-train.hat2sigma2)/true.sigma2]
sigmas[,delta.ape:=abs(true.sigma2-train.deltasigma2)/true.sigma2]

summary(sigmas)

mean.sigmas <- sigmas[,list(true.sigma2=mean(true.sigma2),train.sigma2=mean(train.sigma2),
	train.estsigma2=mean(train.estsigma2),train.adjsigma2=mean(train.adjsigma2),
	train.hatsigma2=mean(train.hatsigma2),train.hat2sigma2=mean(train.hat2sigma2),
	train.deltasigma2=mean(train.deltasigma2)),by=k]

mean.sigmas[,sample.ape:=abs(true.sigma2-train.sigma2)/true.sigma2]
mean.sigmas[,unadj.ape:=abs(true.sigma2-train.estsigma2)/true.sigma2]
mean.sigmas[,adj.ape:=abs(true.sigma2-train.adjsigma2)/true.sigma2]
mean.sigmas[,hat.ape:=abs(true.sigma2-train.hatsigma2)/true.sigma2]
mean.sigmas[,hat2.ape:=abs(true.sigma2-train.hat2sigma2)/true.sigma2]
mean.sigmas[,delta.ape:=abs(true.sigma2-train.deltasigma2)/true.sigma2]

summary(mean.sigmas)

mapes <- lapply(sims,`[[`,"mapes")
mapes <- rbindlist(mapes)
summary(mapes)

r2s <- lapply(sims,`[[`,"r2s")
r2s <- rbindlist(r2s)
r2s[,train.ape:=abs((train-test)/test)]
r2s[,nonstoch.ape:=abs((train-nonstoch)/nonstoch)]
r2s[,nonstoch.hetero.ape:=abs((nonstoch.hetero-nonstoch)/nonstoch)]
r2s[,stoch.nonstoch.ape:=abs((stoch.nonstoch-nonstoch)/nonstoch)]
r2s[,oos.ape:=abs((oos-test)/test)]
saveRDS(r2s,'r2s.homo45_125.RDS')

mean.r2s <- r2s[,list(train=mean(train),test=mean(test),nonstoch=mean(nonstoch),
	nonstoch.hetero=mean(nonstoch.hetero),stoch.nonstoch=mean(stoch.nonstoch),
	oos=mean(oos)),by=k]

mean.r2s[,train.ape:=abs((train-test)/test)]
mean.r2s[,nonstoch.ape:=abs((train-nonstoch)/nonstoch)]
mean.r2s[,nonstoch.hetero.ape:=abs((nonstoch.hetero-nonstoch)/nonstoch)]
mean.r2s[,stoch.nonstoch.ape:=abs((stoch.nonstoch-nonstoch)/nonstoch)]
mean.r2s[,oos.ape:=abs((oos-test)/test)]

write.table(sigmas,"homo45_125.sigmas.txt",sep='\t',row.names=FALSE)
write.table(mean.sigmas,"homo45_125.meansigmas.txt",sep='\t',row.names=FALSE)
write.table(r2s,"homo45_125.r2s.txt",sep='\t',row.names=FALSE)
write.table(mean.r2s,"homo45_125.meanr2s.txt",sep='\t',row.names=FALSE)
write.table(mapes,"homo45_125.mapes.txt",sep='\t',row.names=FALSE)