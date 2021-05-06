CleanData <- function(phy, data) {
  colnames(data) <- c("Species","MigratoryTrait")
  trait <- as.numeric(data$MigratoryTrait)
  names(trait) <- data$Species
  datatree <- treedata(phy, trait)
  return(datatree)
}

VisualizeData <- function(phy, data, path) {
    plotTree(phy,type="fan", fsize=0.0001, ftype="i")
    cols <- setNames(c("#b2182b", "#2166ac"), c("0", "1"))
    tiplabels(pie=as.matrix(data)[,1], piecol=cols,cex=0.2)
    legend("topleft",levels(data),pch=21,pt.bg= cols , pt.cex=2.2)
}

clean.df <- function(data, trait.name){
  data.df <- as.data.frame(data)
  data.df$Species <- rownames(data.df)
  rownames(data.df) <- NULL
  colnames(data.df) <- c(trait.name,"Species")
  data.df[,trait.name] <- as.factor(data.df[, trait.name])
  data.df <- data.df[,c("Species", trait.name)]
  return(data.df)
}
