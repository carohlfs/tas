options(digits=22,scipen=999,stringAsFactors=FALSE)

library(data.table)
library(lmtest)
library(parallel)
library(tidyr)
library(ggplot2)

# Figure 1a: Simulation (Hetero, 45, 80) MSE Comparison

k <- 45
plist <- 2:(k+1)
iters <- 1000

sims <- readRDS('sims.heterosked.RDS')
rss <- rbindlist(rss)

mean.rss <- rss[,list(train=mean(train),test=mean(test),nonstoch=mean(nonstoch),
	nonstoch.hetero=mean(nonstoch.hetero),stoch.nonstoch=mean(stoch.nonstoch),
	oos=mean(oos),loo=mean(loo)),by=k]

Err.train <- mean.rss[,train]
Eerr <- mean.rss[,test]
Ehat <- mean.rss[,oos]
nonstoch <- mean.rss[,nonstoch]
nonstoch.hat <- mean.rss[,nonstoch.hetero]
loo <- mean.rss[,loo]

png(file="Figure_1a.png",width=500,height=400)
par(mar=c(5,4,1,1))
ylims=c(0,1.2)
matplot(plist,Err.train,type="n",col="black",xlab="Model Complexity (df)",ylab="Prediction Error",ylim=ylims,lty=1)
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "gray90") # Color
lines(plist,Err.train,col="black",lwd=4,lty=1)
lines(plist,nonstoch,col="white",lwd=4,lty=4)
lines(plist,nonstoch.hat,col="gray50",lwd=4,lty=2)
lines(plist,Eerr,col="white",lwd=4,lty=1)
lines(plist,Ehat,col="gray50",lwd=4,lty=1)
lines(plist,loo,col="black",lwd=4,lty=4)
legend(x = "topleft",          # Position
        ncol = 2,
       bg = "gray95",
       legend = c("Training","Test (Non-Stochastic)","Proposed (Non-Stochastic)","Test (Stochastic)","Proposed (Stochastic)","PRESS"),  # Legend texts
       lty = c(1, 4, 2, 1, 1, 4),           # Line types
       col = c('black', 'white', 'gray75','white', 'gray50', 'black'),           # Line colors
       lwd = c(4,4,4,4,4,4))                 # Line width
dev.off()

# Figure 1b: Empirical (50%) R-squared Comparison

k <- 50
plist <- 2:(k+1)
iters <- 1000

rss <- readRDS('rss50.RDS')

mean.rss <- rss[,list(train=mean(train),test=mean(test),oos=mean(oos),loo=mean(loo)),by=k]

Err.train <- mean.rss[,train]
Eerr <- mean.rss[,test]
Ehat <- mean.rss[,oos]
loo <- mean.rss[,loo]

png(file="Figure_1b.png",width=500,height=400)
par(mar=c(5,4,1,1))
ylims=c(0,0.65)
matplot(plist,Err.train,type="n",col="black",xlab="Model Complexity (df)",ylab="Prediction Error",ylim=ylims,lty=1)
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "gray90") # Color
lines(plist,Err.train,col="black",lwd=4,lty=1)
lines(plist,Eerr,col="white",lwd=4,lty=1)
lines(plist,Ehat,col="gray50",lwd=4,lty=1)
lines(plist,loo,col="black",lwd=4,lty=4)
legend(x = "topleft",          # Position
        ncol = 2,
       bg = "gray95",
       legend = c("Training","Test","Proposed","PRESS"),  # Legend texts
       lty = c(1, 1, 1, 4),           # Line types
       col = c('black', 'white', 'gray50', 'black'),           # Line colors
       lwd = c(4,4,4,4))                 # Line width
dev.off()

# Figure 2a: Simulation (Hetero 45) Proposed MAPE

sizes45 <- c(50,80,100,125,250,500,1000)
sizes950 <- c(1000,1250,1500,1750,2000,2500,5000)
type <- c("homo","hetero")

files45 <- c(paste0("simulation.data.homo45_",sizes45,".RDS"),paste0("simulation.data.hetero45_",sizes45,".RDS"))
files950 <- c(paste0("simulation.data.homo",sizes950,".RDS"),paste0("simulation.data.hetero",sizes950,".RDS"))

files <- c(files45,files950)

sizes <- gsub(".*_","",files)
sizes <- gsub(".*o","",sizes)
sizes <- gsub(".RDS","",sizes)

types <- gsub("simulation.data.","",files)
types <- gsub("_.*","",types)
types <- gsub(".RDS","",types)
types <- gsub("o(1|2|5).*","o",types)

outputs <- paste0("mapes.",types,"_",sizes,".RDS")

hetero45.mapes <- outputs[grepl("hetero45",outputs)]
hetero45.data <- lapply(hetero45.mapes,readRDS)
mapes <- rbindlist(hetero45.data)
mapes[size=="1000",size:="1,000"]
mapes[,size:=factor(size,levels=c("50","80","100","125","250","500","1,000"))]
setkey(mapes,size,j)

press.mapes <- mapes[,list(size,j+1,press)]
proposed.mapes <- mapes[,list(size,j+1,proposed)]

mape.names <- c("Training Sample Size","Model Complexity (df)","MAPE")
setnames(press.mapes,mape.names)
setnames(proposed.mapes,mape.names)

g <- ggplot(proposed.mapes, aes(x = `Model Complexity (df)`, y = `Training Sample Size`, fill = MAPE)) +
  geom_tile() + scale_fill_continuous(labels = scales::percent,low="white",high="black",limits=c(0,1.75)) + 
  theme(text = element_text(size=40), legend.key.height = unit(2, 'cm'), panel.background = element_rect(fill="white",color="white"),
  	panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "black")) 
png(file="Figure_2a.png",width=800,height=500)
print(g)
dev.off()

# Figure 2b: Simulation (Hetero 45) PRESS MAPE

g <- ggplot(press.mapes, aes(x = `Model Complexity (df)`, y = `Training Sample Size`, fill = MAPE)) +
  geom_tile() + scale_fill_continuous(labels = scales::percent,low="white",high="black",limits=c(0,1.75)) + 
  theme(text = element_text(size=40), legend.key.height = unit(2, 'cm'), panel.background = element_rect(fill="white",color="white"),
  	panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "black")) 
png(file="Figure_2b.png",width=800,height=500)
print(g)
dev.off()

# Figure 2c: Empirical Proposed MAPE

mapes25 <- readRDS('mri25.mapes.RDS')
mapes25[,size:="25%"]
mapes50 <- readRDS('mri50.mapes.RDS')
mapes50[,size:="50%"]
mapes75 <- readRDS('mri75.mapes.RDS')
mapes75[,size:="75%"]
mapes <- rbind(mapes25,mapes50,mapes75)

press.mapes <- mapes[,list(size,j+1,press)]
proposed.mapes <- mapes[,list(size,j+1,proposed)]

mape.names <- c("Training Sample %","Model Complexity (df)","MAPE")
setnames(press.mapes,mape.names)
setnames(proposed.mapes,mape.names)

g <- ggplot(proposed.mapes, aes(x = `Model Complexity (df)`, y = `Training Sample %`, fill = MAPE)) +
  geom_tile() + scale_fill_continuous(labels = scales::percent,low="white",high="black",limits=c(0,1.0)) + 
  theme(text = element_text(size=40), legend.key.height = unit(2, 'cm'), panel.background = element_rect(fill="white",color="white"),
  	panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "black")) 
png(file="Figure_2c.png",width=800,height=500)
print(g)
dev.off()

# Figure 2d: Empirical PRESS MAPE

g <- ggplot(press.mapes, aes(x = `Model Complexity (df)`, y = `Training Sample %`, fill = MAPE)) +
  geom_tile() + scale_fill_continuous(labels = scales::percent,low="white",high="black",limits=c(0,1.0)) + 
  theme(text = element_text(size=40), legend.key.height = unit(2, 'cm'), panel.background = element_rect(fill="white",color="white"),
  	panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "black")) 
png(file="Figure_2d.png",width=800,height=500)
print(g)
dev.off()

# Figure 3a: Simulation (N=80) 2D Error Density
hetero.data <- readRDS('simulation.data.hetero.RDS')

hetero.data[,e2prop:=r+h]
hetero.data <- hetero.data[source=="test",list(e2est,e2prop)]
hetero.data[,e2est:=floor(e2est*2500)/2500]
hetero.data[,e2prop:=floor(e2prop*500)/500]
hetero.data[,count:=as.numeric(1)]
hetero.data <- hetero.data[,list(count = sum(count)),by=list(e2est,e2prop)]
hetero.data[,count:=count*2500*500/450000000]
setkey(hetero.data,e2est,e2prop)

options(digits=1,scipen=999,stringAsFactors=FALSE)

setnames(hetero.data,"e2est","Squared Test Error (Estimated)")
setnames(hetero.data,"e2prop","Squared Test Error (Proposed)")
setnames(hetero.data,"count","Density")

g <- ggplot(hetero.data[`Squared Test Error (Estimated)`<=0.1 & `Squared Test Error (Proposed)`>= -2 & `Squared Test Error (Proposed)`<= 3], aes(x = `Squared Test Error (Proposed)`, y = `Squared Test Error (Estimated)`, fill = Density)) +
  geom_tile() + scale_fill_gradient(low="white",high="black") + theme(text = element_text(size=120), legend.key.height = unit(3, 'cm'),
  	panel.background = element_rect(fill="white",color="white"),panel.grid.major = element_blank())
png(file="Figure_3a.png",width=3000,height=2400)
print(g)
dev.off()

# Figure 3b: Empirical (50%) 2D Error Density
components50 <- readRDS('components50.RDS')
components50 <- lapply(components50,`[[`,"components")
components50 <- rbindlist(components50)

components50[,e2prop:=r+h]
components50 <- components50[source=="test",list(e2est,e2prop)]

components50[,e2est:=floor(e2est*2500)/2500]
components50[,e2prop:=floor(e2prop*500)/500]
components50[,count:=as.numeric(1)]
components50 <- components50[,list(count = sum(count)),by=list(e2est,e2prop)]
components50[,count:=count*2500*500/10627300]
setkey(components50,e2est,e2prop)

options(digits=1,scipen=999,stringAsFactors=FALSE)

setnames(components50,"e2est","Squared Test Error (Estimated)")
setnames(components50,"e2prop","Squared Test Error (Proposed)")
setnames(components50,"count","Density")

g <- ggplot(components50[`Squared Test Error (Estimated)`<=0.1 & `Squared Test Error (Proposed)`>= -2 & `Squared Test Error (Proposed)`<= 3], aes(x = `Squared Test Error (Proposed)`, y = `Squared Test Error (Estimated)`, fill = Density)) +
  geom_tile() + scale_fill_gradient(low="white",high="black") + theme(text = element_text(size=120), legend.key.height = unit(3, 'cm'),
    panel.background = element_rect(fill="white",color="white"),panel.grid.major = element_blank()) + xlim(-2,3)
png(file="Figure_3b.png",width=3000,height=2400)
print(g)
dev.off()

# Figure 4a: Simulation (N=80) Hats

scale.vec <- function(vector){
	return(1000*vector/sum(vector))
}

hetero.data <- readRDS('simulation.data.hetero.RDS')

train.hetero.5 <- hetero.data[source=="training" & j==5,hat]
train.hetero.25 <- hetero.data[source=="training" & j==25,hat]
train.hetero.45 <- hetero.data[source=="training" & j==45,hat]

test.hetero.5 <- hetero.data[source=="test" & j==5,hat]
test.hetero.25 <- hetero.data[source=="test" & j==25,hat]
test.hetero.45 <- hetero.data[source=="test" & j==45,hat]

vectors <- list(train.hetero.5,train.hetero.25,train.hetero.45,test.hetero.5,test.hetero.25,test.hetero.45)
densities <- lapply(vectors,density,bw="SJ",n=65536,from=0,to=65.536)

support <- seq(0.001,65.536,0.001)

dens.graphs <- lapply(densities,`[[`,"y")
dens.graphs <- Reduce(cbind,dens.graphs)
dens.graphs <- data.table(dens.graphs)
density.names <- c("train.hetero.5","train.hetero.25","train.hetero.45","test.hetero.5","test.hetero.25","test.hetero.45")
setnames(dens.graphs,density.names)
dens.graphs[,(density.names):=lapply(.SD,scale.vec),.SDcols=density.names]

png(file="Figure_4a.png",width=500,height=400)
par(mar=c(5,5,1,1))
ylims=c(0,15)
matplot(support[1:2500],dens.graphs[1:2500,train.hetero.5],type="n",col="black",xlab="Leverage",ylab="Density",ylim=ylims,lty=1,
	cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "gray90") # Color
lines(support[1:2500],dens.graphs[1:2500,train.hetero.5],col="black",lwd=5,lty=1)
lines(support[1:2500],dens.graphs[1:2500,train.hetero.25],col="black",lwd=3,lty=3)
lines(support[1:2500],dens.graphs[1:2500,train.hetero.45],col="black",lwd=1,lty=5)
lines(support[1:2500],dens.graphs[1:2500,test.hetero.5],col="white",lwd=5,lty=1)
lines(support[1:2500],dens.graphs[1:2500,test.hetero.25],col="white",lwd=3,lty=3)
lines(support[1:2500],dens.graphs[1:2500,test.hetero.45],col="white",lwd=1,lty=5)
legend(x = "topright",          # Position
        ncol = 2,
       bg = "gray95",
       legend = c("Training, k=5","Training, k=25","Training, k=45","Test, k=5","Test, k=25","Test, k=45"),  # Legend texts
       lty = c(1, 3, 5, 1, 3, 5),           # Line types
       col = c('black','black','black','white','white','white'),           # Line colors
       lwd = c(5,3,1,5,3,1),
       cex=1.5)                 # Line width
dev.off()

# Figure 4b: Empirical (50%) Hats
components50 <- readRDS('components50.RDS')
components50 <- lapply(components50,`[[`,"components")
components50 <- rbindlist(components50)

train50.5 <- components50[source=="training" & j==5,hat]
train50.25 <- components50[source=="training" & j==25,hat]
train50.45 <- components50[source=="training" & j==45,hat]

test50.5 <- components50[source=="test" & j==5,hat]
test50.25 <- components50[source=="test" & j==25,hat]
test50.45 <- components50[source=="test" & j==45,hat]

vectors <- list(train50.5,train50.25,train50.45,test50.5,test50.25,test50.45)
densities <- lapply(vectors,density,bw="SJ",n=65536,from=0,to=65.536)

support <- seq(0.001,65.536,0.001)

dens.graphs <- lapply(densities,`[[`,"y")
dens.graphs <- Reduce(cbind,dens.graphs)
dens.graphs <- data.table(dens.graphs)
density.names <- c("train50.5","train50.25","train50.45","test50.5","test50.25","test50.45")
setnames(dens.graphs,density.names)
dens.graphs[,(density.names):=lapply(.SD,scale.vec),.SDcols=density.names]

png(file="Figure_4b.png",width=500,height=400)
par(mar=c(5,5,1,1))
ylims=c(0,30)
matplot(support[1:1000],dens.graphs[1:1000,train50.5],type="n",col="black",xlab="Leverage",ylab="Density",ylim=ylims,lty=1,
	cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)	
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "gray90") # Color
lines(support[1:1000],dens.graphs[1:1000,train50.5],col="black",lwd=5,lty=1)
lines(support[1:1000],dens.graphs[1:1000,train50.25],col="black",lwd=3,lty=3)
lines(support[1:1000],dens.graphs[1:1000,train50.45],col="black",lwd=1,lty=5)
lines(support[1:1000],dens.graphs[1:1000,test50.5],col="white",lwd=5,lty=1)
lines(support[1:1000],dens.graphs[1:1000,test50.25],col="white",lwd=3,lty=3)
lines(support[1:1000],dens.graphs[1:1000,test50.45],col="white",lwd=1,lty=5)
legend(x = "topright",          # Position
        ncol = 2,
       bg = "gray95",
       legend = c("Training, k=5","Training, k=25","Training, k=45","Test, k=5","Test, k=25","Test, k=45"),  # Legend texts
       lty = c(1, 3, 5, 1, 3, 5),           # Line types
       col = c('black','black','black','white','white','white'),           # Line colors
       lwd = c(5,3,1,5,3,1),
       cex=1.5)                 # Line width
dev.off()

# Figure 5a: Simulation Errors by Leverage (Hetero, 45, 80)
leverage <- c("0.0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0",">1.0")
leverage <- factor(leverage,levels=leverage)

training <- c(0.6238713,0.2981982,0.177047,0.1230406,0.06501982,NA)
press <- c(0.7626606,0.6020672,0.6816414,0.9761453,1.847963,NA)
test <- c(0.7840133,0.6276797,0.5906652,0.608969,0.6562306,0.8095599)
proposed <- c(0.7956379,0.6234469,0.5799454,0.5946818,0.6364819,0.7714633)

png(file="Figure_5a.png",width=500,height=400)
par(mar=c(5,5,1,1))
ylims=c(0,2.0)
matplot(leverage,training,type="n",col="black",xlab="Leverage",ylab="Prediction Error",ylim=ylims,lty=1,,
	cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "gray90") # Color
lines(leverage,training,col="black",lwd=4,lty=1)
lines(leverage,test,col="white",lwd=4,lty=1)
lines(leverage,proposed,col="gray50",lwd=4,lty=1)
lines(leverage,press,col="black",lwd=4,lty=4)
legend(x = "topleft",          # Position
        ncol = 2,
       bg = "gray95",
       legend = c("Training","Test","Proposed",expression(italic(PRESS/n))),  # Legend texts
       lty = c(1, 1, 1, 4),           # Line types
       col = c('black', 'white', 'gray50', 'black'),           # Line colors
       lwd = c(4,4,4,4),
       cex=1.3)                 # Line width
dev.off()

# Figure 5b: Empirical Errors by Leverage (50%)
leverage <- c("0.0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0",">1.0")
leverage <- factor(leverage,levels=leverage)

training <- c(0.3162332,0.2426412,0.2085336,0.1050266,0.07383992,NA)
press <- c(0.3905177,0.4670012,0.7750694,1.049073,4.779605,NA)
test <- c(0.3855432,0.4221322,0.509984,0.6680022,0.7867661,1.289873)
proposed <- c(0.3878625,0.4194544,0.520688,0.6416003,0.7683157,1.307818)

png(file="Figure_5b.png",width=500,height=400)
par(mar=c(5,5,1,1))
ylims=c(0,5.0)
matplot(leverage,training,type="n",col="black",xlab="Leverage",ylab="Prediction Error",ylim=ylims,lty=1,,
	cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "gray90") # Color
lines(leverage,training,col="black",lwd=4,lty=1)
lines(leverage,test,col="white",lwd=4,lty=1)
lines(leverage,proposed,col="gray50",lwd=4,lty=1)
lines(leverage,press,col="black",lwd=4,lty=4)
legend(x = "topleft",          # Position
        ncol = 2,
       bg = "gray95",
       legend = c("Training","Test","Proposed",expression(italic(PRESS/n))),  # Legend texts
       lty = c(1, 1, 1, 4),           # Line types
       col = c('black', 'white', 'gray50', 'black'),           # Line colors
       lwd = c(4,4,4,4),
       cex=1.3)                 # Line width
dev.off()


