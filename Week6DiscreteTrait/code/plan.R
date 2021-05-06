plan <- drake_plan(
  tree = fishtree_phylogeny(type = "chronogram"),
  discrete.data = read.csv("data/MigratoryStatusAll.csv", stringsAsFactors=FALSE), #Freshwater fish migratory species of the amazon
  data.clean =  CleanData(tree, discrete.data),
  VisualizeData(data.clean$phy, data.clean$data, "results/traitmap"),
  cleaned.discrete.phyDat = phangorn::phyDat(data.clean$data, type="USER",levels=c(0,1)), 
  anc.p =  phangorn::ancestral.pars(data.clean$phy, cleaned.discrete.phyDat), 
  phangorn::plotAnc(data.clean$phy, anc.p, pos=1, cex=0.1, cex.pie =0.3, type="fan"),
  anc.ml =  ancestral.pml(pml(data.clean$phy, cleaned.discrete.phyDat), type="ml"),
  phangorn::plotAnc(data.clean$phy, anc.ml, pos=1, cex=0.1, cex.pie =0.3, type="fan"),
  
  ########### BIOLOGICAL QUESTIONS ##################
  #How does this differ from parsimony?#
  
  # Parsimony provide absolute/discrete reconstruction of ancestral states, while maximum likelihood shows support 
  #in terms of probability represented in pie charts
  
  ##Why does it differ from parsimony?##
  
  # Parsimony reconstruction deals with
  
  ##What does uncertainty mean?##
  
  #Uncertainty is reflected as the degree of support for the reconstruction of each one of the nodes, which,
  #is impacted from the limitation in reconstructions in more ancestral nodes. Uncertainty decreases as we go deeper in the tree.
  
  ##1. How can you estimate transition rates between states? Do it.##
  
  CalculateTransitionRates =  corHMM::corHMM(data.clean$phy, clean.df(data.clean$data, "Migration"), rate.cat = 1),
  
  #2. How could you examine if transition rates are equal?##
  
  TransitionRates = print(CalculateTransitionRates$solution),
  
  # Transition rates are not equal as there is a higher transition rate from  migratory (1) to non migratory (0),
  # than non-migratory (0) to migratory (1)
  
  # 3. Think about the Lewis (2001) MKV model. Are your traits all variable? Will using this make sense for your data? Try using it. Do results change?###
  
  # Migratory trait is variable as it is present across multiple clades across the tree (no autapomorphies) and both states occur. 
  #For that reason, seems reasonable to use Lewis (2001) MKV model. However, I do not notice major differences visually with the previous approach.
  
  anc.Mkv =  phangorn::ancestral.pml(pml(data.clean$phy, cleaned.discrete.phyDat, Mkv = TRUE), type="ml"), #MKV through Maximum Likelihood
  phangorn::plotAnc(data.clean$phy, anc.Mkv, 1, cex=0.1, size=0.0001, cex.pie =0.3, type="fan")
  #  4. How could you test order of state evolution?
  )
