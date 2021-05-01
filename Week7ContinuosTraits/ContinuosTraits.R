library(ape) #utility fns
library(geiger) #utilty fns
library(OUwie)
library(fishtree)
library(dplyr)
library(tidyr)

continuous.data <- read.csv("D:/OneDrive - Pontificia Universidad Javeriana/Documents/NutritionalValuesFishes/Nutrients_EntireDataSet.csv")
tree <- fishtree_phylogeny(type="chronogram")

CleanData <- function(data, phylo) {
 
  data$Species <- gsub(' ', '_' , continuous.data$Species)
  data <- data %>%
    group_by(Species) %>%
    summarize(mean_Pr = mean(Pr)) %>%
    drop_na()  %>%
    filter(Species %in% phylo$tip.label)
  
  trait <- as.numeric(data$mean_Pr)
  names(trait) <- data$Species
  
  datatree <- treedata(tree, trait)
  return(datatree)
 
}

data <- CleanData(continuous.data, tree)

VisualizeData <- function(phy, data) {
  phytools::contMap(phy, data)
}

#### This is how is fitted here
BM1 <- geiger::fitContinuous(data$phy, data$data, model="BM")
print(paste("The rate of evolution is", BM1$opt$sigsq, "in units of", "time"))


OU1 <- fitContinuous(data$phy, data$data, model="OU")
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
reconstruction.info <- ace(one.discrete.char, data$phy, type="discrete", method="ML", CI=TRUE)
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

cont.data <- data.frame(species=rownames(data$data), regime=best.states, trait=data$data[,1])

labeled.tree <- data$phy
labeled.tree$node.label <- best.states
nodeBased.OUMV <- OUwie(data$phy, cont.data, model="OUMV", simmap.tree=FALSE, diagn=FALSE)
print(nodeBased.OUMV)

models <- c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA")
results <- lapply(models, RunSingleOUwieModel, phy=tree, data=trait)

AICc.values<-sapply(results, "[[", "AICc")
names(AICc.values)<-models
AICc.values<-AICc.values-min(AICc.values)

print(AICc.values) #The best model is the one with smallest AICc score

best<-results[[which.min(AICc.values)]] #store for later

print(best) #prints info on best model
