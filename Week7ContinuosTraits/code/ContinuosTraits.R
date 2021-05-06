source("code/functions.R")
source("code/packages.R")


#### CLEANING DATA ####
continuous.data <- read.csv("D:/OneDrive - Pontificia Universidad Javeriana/Documents/NutritionalValuesFishes/data/Nutrients_EntireDataSet_AllSpecies.csv")
tree <- fishtree_phylogeny(type="chronogram")
data <- CleanData(continuous.data, tree)

#### VISUALIZE DATA ####

VisualizeData(data$phy, data$data)
#### This is how is fitted here
BM1 <- geiger::fitContinuous(data$phy, data$data, model="BM")
print(paste("The rate of evolution is", BM1$opt$sigsq, "in units of", "time"))

OU1 <- fitContinuous(phytools::force.ultrametric(data$phy), data$data, model="OU")
par(mfcol(c(1,2)))
plot(tree, show.tip.label=FALSE)
ou.tree <- rescale(tree, model="OU", OU1$opt$alpha)
plot(ou.tree)

AIC.BM1 <- AIC(BM1)
AIC.OU1 <- AIC(OU1)
delta.AIC.BM1 <- AIC.BM1 - min(AIC.BM1, AIC.OU1) 
delta.AIC.OU1 <- AIC.OU1 - min(AIC.BM1, AIC.OU1) 

##OUwie runs##

one.discrete.char <- as.numeric(data$data < mean(data$data))
reconstruction.info <- ace(one.discrete.char, phytools::force.ultrametric(data$phy), type="discrete", method="ML", CI=TRUE)
best.states <- colnames(reconstruction.info$lik.anc)[apply(reconstruction.info$lik.anc, 1, which.max)]

clean.df <- function(data, trait.name){
  data.df <- as.data.frame(data)
  data.df$Species <- rownames(data.df)
  rownames(data.df) <- NULL
  colnames(data.df) <- c(trait.name,"Species")
  data.df[,trait.name] <- as.factor(data.df[, trait.name])
  data.df <- data.df[,c("Species", trait.name)]
  return(data.df)
}

cont.data <- data.frame(species=rownames(data$data), regime=one.discrete.char, trait=data$data[,1])

labeled.tree <- data$phy
labeled.tree$node.label <- best.states

nodeBased.OUMV <- OUwie(labeled.tree, cont.data, model="OUMV", simmap.tree=FALSE, diagn=FALSE)
print(nodeBased.OUMV)
# What do the numbers mean?

RunSingleOUwieModel<-function(phy, data, models){
  mod<-OUwie(phy, data, models, simmap.tree = F, diagn = F, root.age = NULL)
}

models <- c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA")
results <- lapply(models, RunSingleOUwieModel, phy=labeled.tree, data=cont.data)

AICc.values<-sapply(results, "[[", "AICc")
names(AICc.values)<-models
AICc.values<-AICc.values-min(AICc.values)

print(AICc.values) #The best model is the one with smallest AICc score

best<-results[[which.min(AICc.values)]] #store for later

print(best) #prints info on best model


alpha.values<-seq(from= 0.0000000000001 , to= 0.01 , length.out=50)

likelihood.values <- rep(NA, length(alpha.values))
for (iteration in sequence(length(alpha.values))) {
  likelihood.values[iteration] <- OUwie.fixed(labeled.tree, cont.data, model="OUMV", alpha=rep(alpha.values[iteration],2), sigma.sq=best$solution[2,], theta=best$theta[,1])$loglik 
}

plot(x= alpha.values, y= likelihood.values, xlab="alpha", ylab="likehlihood", type="l", bty="n", ylim = c(2*best$loglik, best$loglik))
points(x=best$solution[1,1], y=best$loglik, pch=16, col="red")
text(x=best$solution[1,1], y=best$loglik, "unconstrained best", pos=4, col="red")
abline(h=2*best$loglik, lty="dotted") #Two log-likelihood

#Now, let’s try looking at both theta parameters at once, keeping the other parameters at their MLEs

nreps<-400
theta1.points<-c(best$theta[1,1], rnorm(nreps-1, best$theta[1,1], 5*best$theta[1,2])) #center on optimal value, have extra variance
theta2.points<-c(best$theta[2,1], rnorm(nreps-1, best$theta[2,1], 5*best$theta[2,2])) #center on optimal value, have extra variance
likelihood.values<-rep(NA,nreps)

for (iteration in sequence(nreps)) {
  likelihood.values[iteration] <- OUwie.fixed(labeled.tree, cont.data, model="OUMV", alpha=best$solution[1,], sigma.sq=best$solution[2,], theta=c(theta1.points[iteration], theta2.points[iteration]))$loglik
}

likelihood.differences<-(-(likelihood.values-max(likelihood.values)))

interpolated.points<-interp(x=theta1.points, y=theta2.points, z= likelihood.differences, linear=FALSE, extrap=TRUE, xo=seq(min(theta1.points), max(theta1.points), length = 400), yo=seq(min(theta2.points), max(theta2.points), length = 400))

contour(interpolated.points, xlim=range(c(theta1.points, theta2.points)),ylim=range(c(theta1.points, theta2.points)), xlab="Theta 1", ylab="Theta 2", levels=c(2,5,10),add=FALSE,lwd=1, bty="n", asp=1)

points(x=best$theta[1,1], y=best$theta[2,1], col="red", pch=16)


points(x=cont.data$trait[which(cont.data$regime==0)],y=rep(min(c(theta1.points, theta2.points)), length(which(cont.data$regime==0))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 1, plotted along x axis
points(y=cont.data$trait[which(cont.data$regime==1)],x=rep(min(c(theta1.points, theta2.points)), length(which(cont.data$regime==1))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 2, plotted along y axis

trait.ordered<-data.frame(cont.data[,2], cont.data[,2],row.names=cont.data[,1])
trait.ordered<- trait.ordered[labeled.tree$tip.label,]
z<-trait.ordered[,1]
names(z)<-rownames(trait.ordered)
tree.mapped<-make.simmap(labeled.tree,z,model="ER",nsim=1)
leg<-c("black","red")
names(leg)<-c(1,2)
plotSimmap(tree.mapped,leg,pts=FALSE,ftype="off", lwd=1)

simmapBased<-OUwie(tree.mapped,cont.data,model="OUMV", simmap.tree=TRUE, diagn=FALSE)
print(simmapBased)
print(best)
