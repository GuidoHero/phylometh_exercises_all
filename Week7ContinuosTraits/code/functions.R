
RunSingleOUwieModel<-function(phy, data, models){
  mod<-OUwie(phy, data, models, simmap.tree = F, diagn = F, root.age = NULL)
}


CleanData <- function(data, phylo) {
  data$Species <- gsub(' ', '_' , data$Species)
  
  data <- data %>% 
    filter(Food.Type == "FreshFish") %>% 
    filter(Tissue== "Edible parts")  %>% 
    group_by(Species) %>%
    summarize(mean_Pr = mean(Pr, na.rm=T)) %>%
    drop_na()  %>%
    filter(Species %in% phylo$tip.label) 
  
  trait <- as.numeric(data$mean_Pr)
  names(trait) <- data$Species
  
  datatree <- geiger::treedata(tree, as.data.frame(trait))
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