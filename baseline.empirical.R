require(data.table)
library(lmtest)
require(parallel)
require(tidyr)
# library(latex2exp)

cores <- 20

set.seed(04251978)

formula.j <- function(j){
    paste("y",paste0("V",1:j,collapse=' + '),sep=' ~ ')
}

train.data <- function(j,dataset){
    eqn <- formula.j(j)
    reg <- lm(eqn,data=dataset)
    return(reg)
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

calc.rss <- function(reg,dataset){
    yhat <- as.numeric(predict(reg,dataset))
    ybar <- mean(dataset[,y])
    e <- dataset[,y]-yhat
    rss <- mean(e^2)
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

    yhat <- as.numeric(predict(reg,dataset))
    ybar <- mean(dataset[,y])
    e2est <- (dataset[,y]-yhat)^2
    tss.est <- (dataset[,y]-ybar)^2

    components <- data.table(j=reg$rank-1,hat=hO,e2est=e2est,tss.est=tss.est,ehat2=sqs,r=sigma2hat,h=h2sigma2hat)
    return(components)
}


run.random <- function(i,pct){

    cat(paste0("iteration #",i,"\n"))
    pca[,random:=runif(nrow(pca))]
    train.set <- pca[random<=pct]
    test.set <- pca[random>pct]

    k.regs <- lapply(factors,train.data,train.set)

    train.rss <- sapply(k.regs,calc.rss,train.set)
    test.rss <- sapply(k.regs,calc.rss,test.set)
    oos.rss <- unlist(lapply(k.regs,calc.stoch.rss,test.set))
    loo.rss <- sapply(k.regs,calc.loo.rss)
    rss <- data.table(i=i,k=factors,train=train.rss,test=test.rss,oos=oos.rss,loo=loo.rss)
    
    test.components <- lapply(k.regs,calc.components,test.set)
    test.components <- rbindlist(test.components)
    test.components[,source:="test"]

    train.components <- lapply(k.regs,calc.components,train.set)
    train.components <- rbindlist(train.components)
    train.components[,source:="training"]

    components <- rbind(train.components,test.components)
    components[,obs:=.I]
    components[,rep:=i]

    run.data <- list("rss"=rss,"components"=components)

    return(run.data)
}

pca <- readRDS('mri.RDS')

iters <- 1000

factors <- seq(1,50,1)
plist <- factors+1

pct <- 0.50
runs <- mclapply(1:iters,run.random,pct,mc.cores=cores)
saveRDS(runs,"components50.RDS")

rss <- mclapply(runs,`[[`,"rss")
rss <- rbindlist(rss)
saveRDS(rss,'rss50.RDS')

train <- rss[,list(i,k,train)]
train <- spread(train,k,train)
train[,i:=NULL]
train <- as.matrix(train)
err <- train

test <- rss[,list(i,k,test)]
test <- spread(test,k,test)
test[,i:=NULL]
test <- as.matrix(test)
Err.train <- test

trainoos <- rss[,list(i,k,oos)]
trainoos <- spread(trainoos,k,oos)
trainoos[,i:=NULL]
trainoos <- as.matrix(trainoos)
Err.hat <- trainoos

png(file="MRI_50pct_rss.png",width=500,height=400)
par(mar=c(5,4,1,1))
ylims=c(0,2)
Eerr=apply(err,2,mean)
Ehat=apply(Err.hat,2,mean)
matplot(plist,t(Err.train),type="n",col="white",xlab="Model Complexity (df)",ylab="Prediction Error",ylim=ylims,lty=1)
col=rep(c("black","white","gray50"),c(iters,iters,iters))
bigmat=cbind(t(err),t(Err.train),t(Err.hat))
o=sample(seq(3*iters))
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "gray90") # Color
matlines(plist,bigmat[,o],col=col[o],lty=1,lwd=.5)
lines(plist,apply(Err.train,2,mean),col="black",lwd=4,lty=1)
lines(plist,Eerr,col="gray50",lwd=4,lty=1)
lines(plist,Ehat,col="white",lwd=4,lty=1)
legend(x = "topleft",          # Position
        ncol = 2,
       bg = "gray95",
       legend = c("Test Error","Est. Test Error","Training Error",
       "Mean Test Error","Mean Est. Test Error","Mean Training Error"),  # Legend texts
       lty = c(1, 1, 1, 1, 1, 1),           # Line types
       col = c('white','gray50','black','black', 'white', 'gray50'),           # Line colors
       lwd = c(0.75,0.75,0.75,4,4,4))                 # Line width
dev.off()

rss[,train.ape:=abs((train-test)/test)]
rss[,oos.ape:=abs((oos-test)/test)]
rss[,loo.ape:=abs((loo-test)/test)]

median.rss <- rss[,list(train=median(train),test=median(test),oos=median(oos),loo=median(loo)),by=k]
mean.rss <- rss[,list(train=mean(train),test=mean(test),oos=mean(oos),loo=mean(loo)),by=k]

mean.rss[,train.ape:=abs((train-test)/test)]
mean.rss[,oos.ape:=abs((oos-test)/test)]
mean.rss[,loo.ape:=abs((loo-test)/test)]