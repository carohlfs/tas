
library(data.table)
library(lmtest)
library(tidyr)

n.train <- 80
n.test <- 10000
k <- 45
esigma <- 75
reps <- 1000
alpha <- sqrt((150^2 -75^2)/(sum((1:k)^2)+150^2))
ex.scale <- 1/alpha
y.scale <- sqrt(sum((1:k)^2) + sum(((1:k)/ex.scale)^2) +esigma^2)

set.seed(04251978)

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
	x[,y:=rowSums(.SD),.SDcols=xnames]
	x[,rando:=rnorm(n)]
	x[,e:=rando*(esigma+y/ex.scale)]
	x[,y:=(1/y.scale)*(y+e)]
	return(x)
}

gen.conditional.variance <- function(j,dataset){
	xnames <- paste0("x.",1:j)
	sqs <- paste0("sq.",1:j)
	d <- dataset[,mget(xnames)]
	d[,(sqs):=lapply(.SD,`^`,2),.SDcols=xnames]
	d[,sqsum:=rowSums(.SD),.SDcols=sqs]
	alpha.sqsum <- d[,sqsum]*alpha^2/y.scale^2
	uncondit.portion <- sum(((46-j):1)^2)*(1+alpha^2)/y.scale^2
	return(alpha.sqsum+uncondit.portion)
}

gen.nonstoch <- function(j,train.data,test.data){
	copies <- floor(nrow(test.data)/nrow(train.data))
	used.x <- m.copies.df(train.data[,1:j],copies)
	alt.data <- cbind(used.x,test.data[,(j+1):ncol(test.data)])
	xnames <- paste0("x.",1:k)
	alt.data[,y:=rowSums(.SD),.SDcols=xnames]
	alt.data[,e:=rando*(esigma+y/ex.scale)]
	alt.data[,y:=(1/y.scale)*(y+e)]
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
	sigma2 <- e^2/(1-hat)
	dataset[,y:=sigma2]
	eqn <- formula.j(ncol(dataset)-1)
	sigma.aux <- lm(eqn,data=dataset)
	sigma2 <- hat*e^2/(1-hat)
	dataset[,y:=sigma2]
	sigmahat.aux <- lm(eqn,data=dataset)
	return(list(sigma.aux=sigma.aux,sigmahat.aux=sigmahat.aux))
}

gen.true.resid <- function(j,dataset){
	xvars <- names(dataset)[1:j]
	yhat <- dataset[,(1/y.scale)*rowSums(.SD),.SDcols=xvars]
	e <- dataset[,y]-yhat
	return(e)
}

calc.rss <- function(reg,dataset){
	yhat <- as.numeric(predict(reg,dataset))
	ybar <- mean(dataset[,y])
	e <- dataset[,y]-yhat
	rss <- mean(e^2)
	return(rss)
}

calc.nonstoch.rss <- function(reg){
	yhat <- as.numeric(reg$fitted.values)
	ybar <- mean(yhat)
	e <- as.numeric(reg$residuals)
	hat <- as.numeric(lm.influence(reg)$hat)
	rss <- mean((1+hat)*e^2/(1-hat))
	return(rss)
}

calc.stoch.rss <- function(reg,dataset){

  xmatrix <- data.table(reg$model)
  y <- xmatrix[,y]
  xnames <- copy(names(xmatrix))
  xnames[1] <- "const"
  setnames(xmatrix,xnames)
  xmatrix[,const:=1]
  xmatrix <- as.matrix(xmatrix)
	xpxinv <- chol2inv(qr.R(reg$qr))

  xO <- dataset[,xnames[2:length(xnames)],with=FALSE]
  xO[,const:=1]
  setcolorder(xO,xnames)
  xO <- as.matrix(xO)
  hO <- xO %*% xpxinv
  hO <- hO %*% t(xmatrix)

  hat <- as.numeric(lm.influence(reg)$hat)
  e <- as.numeric(reg$residuals)

  sigma2 <- e^2/(1-hat)
  sigma2 <- as.matrix(sigma2)
  sigma2hat <- as.numeric(hO %*% sigma2)
  h2sigma2hat <- as.numeric(hO^2 %*% sigma2)

  rss <- mean(sigma2hat)
  hss <- mean(h2sigma2hat)

  output <- rss+hss

  return(output)
}


calc.loo.rss <- function(reg){
   y <- reg$model$y
   ybar <- mean(y)
   e <- as.numeric(reg$residuals) 
   hat <- as.numeric(lm.influence(reg)$hat)
   eadj <- e/(1-hat)
   rss <- mean(eadj^2)
   return(rss)
}

calc.components <- function(reg,dataset){
	xmatrix <- data.table(reg$model)
	y <- xmatrix[,y]
	xnames <- copy(names(xmatrix))
	xnames[1] <- "const"
	setnames(xmatrix,xnames)
	xmatrix[,const:=1]
	xmatrix <- as.matrix(xmatrix)
	xpxinv <- chol2inv(qr.R(reg$qr))

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

	hO <- rowSums(hO^2)
	e2 <- gen.true.resid(reg$rank-1,dataset)^2

	yhat <- as.numeric(predict(reg,dataset))
	ybar <- mean(dataset[,y])
	e2est <- (dataset[,y]-yhat)^2
	tss.est <- (dataset[,y]-ybar)^2

	components <- data.table(j=reg$rank-1,e2=e2,hat=hO,e2est=e2est,tss.est=tss.est,ehat2=sqs,r=sigma2hat,h=h2sigma2hat)
	return(components)
}

run.simulation <- function(iteration){
	cat(paste0("rep ",iteration,"\n"))
	training.data <- gen.data(n.train,k,esigma)
	test.data <- gen.data(n.test,k,esigma)
	nonstoch.data <- lapply(1:k,gen.nonstoch,training.data,test.data)

	k.regs <- lapply(1:k,train.data,training.data)

	train.rss <- sapply(k.regs,calc.rss,training.data)
	nonstoch.rss <- mapply(calc.rss,k.regs,nonstoch.data)
	test.rss <- sapply(k.regs,calc.rss,test.data)

	nonstoch.hetero.rss <- sapply(k.regs,calc.nonstoch.rss)
	oos.rss <- sapply(k.regs,calc.stoch.rss,test.data)
	loo.rss <- sapply(k.regs,calc.loo.rss)

	rss <- data.table(i=iteration,k=1:k,train=train.rss,nonstoch=nonstoch.rss,test=test.rss,
		nonstoch.hetero=nonstoch.hetero.rss,oos=oos.rss,loo=loo.rss)

	test.components <- lapply(k.regs,calc.components,test.data)
	test.components <- rbindlist(test.components)
	test.components[,source:="test"]

	train.components <- lapply(k.regs,calc.components,training.data)
	train.components <- rbindlist(train.components)
	train.components[,source:="training"]

	components <- rbind(train.components,test.components)
	components[,obs:=.I]
	components[,rep:=iteration]

	write.table(components,"simulation.data.hetero.txt",sep='\t',row.names=FALSE,col.names= iteration==1,append= iteration!=1)

	return(rss)

}

sims <- lapply(1:reps,run.simulation)
saveRDS(sims,'sims.heterosked.RDS')
rss <- rbindlist(sims)


