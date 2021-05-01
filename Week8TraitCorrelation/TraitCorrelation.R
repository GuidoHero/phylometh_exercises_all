library(geiger)
library(ape)
library(corHMM)

##Continuous data

tree.primates <- read.tree(text="((((Homo:0.21,Pongo:0.21):0.28,Macaca:0.49):0.13,Ateles:0.62):0.38,Galago:1.00);") #using examples from ape ?pic
X <- c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968)
Y <- c(4.74493, 3.33220, 3.36730, 2.89037, 2.30259)
names(X) <- names(Y) <- c("Homo", "Pongo", "Macaca", "Ateles", "Galago")
pic.X <- pic(X, tree.primates)
pic.Y <- pic(Y, tree.primates)

##Now, positivitize the contrasts and do a regression through the origin.

##Discrete data

data(primates)
ls()
print(primates)
require(phytools)

#Just to make this a better dataset for our exercise, I'm going to change some of the states 
#(I want to have all four trait combinations present). For actual analyses, of course, DO NOT INVENT YOUR DATA.

#First, a review of discrete state models:

primates$trait[which(grepl("Hylobates",primates$trait[,1])),2]<-1

trait1<-primates$trait[,2]
names(trait1) <- primates$trait[,1]
primates$tree <- ape::multi2di(primates$tree)
plotSimmap(make.simmap(primates$tree, trait1), pts=FALSE, fsize=0.8)
rate.mat.er <- corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=2, model="ER")
print(rate.mat.er)

#What does this matrix mean?
  # This is a matrix of equal transition rates model for two states in one character

pp.er <- corHMM(primates$tree,primates$trait[,c(1,2)], rate.cat=1, rate.mat=rate.mat.er, node.states="marginal")
print(pp.er)

# What do these results mean?

# This are the transition rates estimates of character states from primate data
# under the assumption of the ER model.

rate.mat.ard<-corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=2, model="ARD")
print(rate.mat.ard)

# And these?
# This is the matrix from an all-rates-different model assuming the the transition rate could be different.

pp.ard<-corHMM:::corHMM(primates$tree,primates$trait[,c(1,2)],rate.cat=1,rate.mat=rate.mat.ard, node.states="marginal")
print(pp.ard)

#which model is better?

# Following AIC, ARD model seems better. However the differences between both are under less than a delta AIC of 2, which can be
# considered a small difference.

pp.er$AIC
pp.ard$AIC
  
#Now let's look at multiple traits.

#This is a matrix with four states

rate.mat.er.4state <-corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=4, model="ER")
print(rate.mat.er.4state)

fourstate.trait<-rep(NA,Ntip(primates$tree))
for(i in sequence(Ntip(primates$tree))) {
  if(primates$trait[i,2]==0 && primates$trait[i,3]==0) {
    fourstate.trait[i]<-0
  }
  if(primates$trait[i,2]==0 && primates$trait[i,3]==1) {
    fourstate.trait[i]<-1
  }
  if(primates$trait[i,2]==1 && primates$trait[i,3]==0) {
    fourstate.trait[i]<-2
  }
  if(primates$trait[i,2]==1 && primates$trait[i,3]==1) {
    fourstate.trait[i]<-3
  }
}

fourstate.data<-data.frame(Genus_sp=primates$trait[,1], T1=fourstate.trait)

print(rayDISC(primates$tree, fourstate.data, ntraits=1, model="ER", node.states="marginal"))
print(rayDISC(primates$tree, fourstate.data, ntraits=1, rate.mat=rate.mat.er.4state, node.states="marginal", model="ARD"))
rate.mat.ard.4state<-corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=4, model="ARD")
print(rate.mat.ard.4state)


#Now let's make the equivalent of a GTR matrix:

rate.mat.gtr.4state<-rate.mat.ard.4state
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(1,4))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(2,6))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(3,8))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(4,6))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(5,7))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(6,7))
print(rate.mat.gtr.4state)

#Now make a model like Pagel 1994

print(rayDISC(primates$tree, fourstate.data, ntraits=1, rate.mat= rate.mat.gtr.4state, node.states="marginal", model="ARD"))

print(corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=2, nstates=2, model="ARD"))
rate.mat.pag94 <- corHMM:::rate.par.drop(rate.mat.ard.4state, drop.par=c(3,5,8,10))

#Now that you have some introduction, there are two routes:

#Construct a model to test if state 1 can never be lost

rate.mat.NoOneLoss <- corHMM:::rate.par.drop(rate.mat.ard, drop.par=1) # Set 2 (1) to 1(0) transition to NA
rate.mat.NoOneLossNoNA <- rate.mat.NoOneLoss
rate.mat.NoOneLossNoNA[is.na(rate.mat.NoOneLossNoNA)] <- 0
primates$trait2 <- primates$trait
primates$trait2[,2] <- primates$trait2[,2] + 1
pp.NoOneLoss <- corHMM(primates$tree,primates$trait2[,c(1,2)], rate.cat=1, rate.mat=rate.mat.NoOneLossNoNA, node.states="none")

#Experiment with the effects of frequencies at the root.

pp.ard.maddfitz <-corHMM:::corHMM(primates$tree,primates$trait[,c(1,2)],rate.cat=1,rate.mat=rate.mat.ard, node.states="marginal", root.p="maddfitz")
pp.ard.yang <-corHMM:::corHMM(primates$tree,primates$trait[,c(1,2)],rate.cat=1,rate.mat=rate.mat.ard, node.states="marginal", root.p="yang")

pp.ard.maddfitz$AIC
pp.ard.yang$AIC 

#Create and use a model to see if transitions from 00 go to 11 only via 01.

# Following AIC, allowing all paths seems better fir model than restricting transition 00-11 via 01. 
# However the differences between both are under less than a delta AIC of 2, which can be
# considered a "non-significant" difference.

rate.mat.ard.pag94.OnlyVia01 <-corHMM:::rate.par.drop(rate.mat.pag94, drop.par=c(5))
rate.mat.ard.All.OnlyVia01 <-corHMM:::rate.par.drop(rate.mat.gtr.4state, drop.par=c(2,3))

pp.All<- rayDISC(primates$tree, fourstate.data, ntraits=1, rate.mat= rate.mat.ard.4state, node.states="marginal", model="ARD")
pp.pag94 <- rayDISC(primates$tree, fourstate.data, ntraits=1, rate.mat= rate.mat.pag94, node.states="marginal", model="ARD")
pp.All.Onlyvia01 <- rayDISC(primates$tree, fourstate.data, ntraits=1, rate.mat= rate.mat.ard.All.OnlyVia01, node.states="marginal", model="ARD")
pp.pag94.Onlyvia01 <- rayDISC(primates$tree, fourstate.data, ntraits=1, rate.mat= rate.mat.ard.pag94.OnlyVia01, node.states="marginal", model="ARD")

print(pp.All$AIC)
print(pp.pag94$AIC)
print(pp.All.Onlyvia01$AIC)
print(pp.pag94.Onlyvia01$AIC)

