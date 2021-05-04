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
  cirplot.trait <- gheatmap(p8, trait, width=.15, colnames = F, color=NULL) + 
                  scale_fill_viridis_c(option="D") 
  return(cirplot.trait)
  }