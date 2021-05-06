CleanData <- function(phy, data) {
  colnames(data) <- c("Species","MigratoryTrait")
  trait <- as.numeric(data$MigratoryTrait)
  names(trait) <- data$Species
  datatree <- treedata(phy, trait)
  return(datatree)
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

VisualizeData <- function (tree, trait){
  circplot <- ggtree(tree, layout = "circular")
  cirplot.trait <- gheatmap(circplot , trait, width=.15, colnames = F, color=NULL) + 
                  scale_fill_viridis_c(option="D") 
  return(plot(cirplot.trait))
  }