#This script extracts global topological network measures from different approaches for every fMRI task type and subject, and collates them together into a concise data.frame

library(PHfconn)
library(igraph)
library(ape)

for (taskFC in c('REST3T','REST7T','MOTOR','GAMBLING','LG','WM','SOCIAL','RELATIONAL','MOVIE'))
{
    
  ######################
  
  #load lists of FC matrices, each individual matrix produced by produced by FCtools package and concatenated in a single cohort matrix: N subject x E edges
  FClist=readRDS(paste0("HCP/taskFC/HCP_",taskFC,"_FC.rds"))

  #get matching participant data (participants ID are stored in the FC matrix)
  beh_data <- read.csv('HCP/HCP_behdata.csv', row.names = 1)
  colnames(beh_data)[1]='participant_id'
  beh_data=beh_data[which(beh_data$participant_id %in% row.names(FClist)),]
  
  #prepare data.frame to store network measures
  NTW_measures=data.frame('GLOB_EFF0.1'=NA, 'CLUSTERING0.1'=NA,
                          'GLOB_EFF0.2'=NA, 'CLUSTERING0.2'=NA,
                          'GLOB_EFF0.3'=NA, 'CLUSTERING0.3'=NA,
                          'PH_BS'=NA,'PH_BD'=NA,'PH_CS'=NA,
                          'MST_LEAF'=NA,'MST_DIAM'=NA,
                          'RAW_FC'=NA)[rep(1, nrow(beh_data)), ]
  row.names(NTW_measures)=beh_data$participant_id
  
  for (sub in 1:nrow(beh_data))
  {
    #extract adjacency matrix from the cohort-wise data
    #matrices were 219x219 as FCtools includes 19 aseg subcortical regions 
    FCmat = matrix(0, nrow = 219, ncol = 219)
    FCmat[upper.tri(FCmat, diag = F)] = as.numeric(FClist[sub,])
    tm <- t(FCmat)
    FCmat[lower.tri(FCmat)] <- tm[lower.tri(tm, diag = F)]
    diag(FCmat) <- 1 #diag needed for PHfconn/MST building, will be ignored by igraph
    FCmat[FCmat<0]=0 #exclude negative edges
    FCmat=FCmat[1:200,1:200] #no aseg
    
    ############################################################
    ####################Basic raw FC############################
    #mean Fisher z‐transformed correlation coefficient of all pairs of nodes
    #fisher z transformation was already applied
    NTW_measures[beh_data$participant_id[sub],'RAW_FC']=mean(FCmat) 
    
    ############################################################
    ####################Graph theory measures###################
    #multiple density thresholds 
    for (thresh in c(.10,.20,.30))
    {
      #apply threshold
      FCmat_t=FCmat
      FCmat_t[FCmat_t<=quantile(FCmat, 1-thresh)]=0
      graphobj=igraph::graph_from_adjacency_matrix(FCmat_t, 
                                                   mode = "undirected", 
                                                   weighted = TRUE, 
                                                   diag=FALSE)
      #extract graph measure
      GLOB_EFF=igraph::global_efficiency(graphobj)
      CLUSTERING=igraph::transitivity(graphobj)
      NTW_measures[beh_data$participant_id[sub],
                   paste0(c('GLOB_EFF','CLUSTERING'),thresh)]=c(GLOB_EFF,CLUSTERING)
    }
    
    ############################################################
    ####################Persistent Homology measures############
    #See definitions from doi: 10.1002/hbm.26304
    #https://pmc.ncbi.nlm.nih.gov/articles/PMC10203816/
    #https://github.com/hyunnamryu/PHfconn/tree/main
    
    #No need to reverse weights so edges get added stronger to weaker in MST building, PHfconn's connectivity_weights_set() function already handles that.
    PHlist=PHfconn::PH_meas(FCmat)[c('BS','BD','CS')]
    NTW_measures[beh_data$participant_id[sub],
                 paste0('PH_',names(PHlist))]=unlist(PHlist) 
    
    ############################################################
    #######################MST measures#########################
    #See definitions from doi: 10.1162/netn_a_00245
    #https://pmc.ncbi.nlm.nih.gov/articles/PMC9207994/#sec15
    
    #kept single-linkage method for all MST-based computation for consistency
    MST=stats::hclust(stats::as.dist(1-FCmat),method="single")
    #recreate matrix with connected nodes of the MST:
    #deduces which edges were selected based on their weight value 
    #(only works if there are no duplicates in the FC matrix)
    uppertrg=FCmat[upper.tri(FCmat)]
    if (length(which(duplicated(uppertrg[which(uppertrg!=0)])==TRUE))!=0)
    {stop("FC matrix has duplicates. hclust can't give edge indices")}
    #Build binary MST adjacency matrix
    MST_matrix <- matrix(0, nrow = nrow(FCmat), ncol = ncol(FCmat))
    #deduce edges from height values
    MST_matrix[which((1 - FCmat) %in% MST$height)] = 1
    MST=MST_matrix
    #Leaf Fraction
    MST_LEAF= sum(rowSums(MST) == 1) / length(rowSums(MST))
    #Diameter (largest distance / total )
    diameter=igraph::diameter(igraph::graph_from_adjacency_matrix(MST, 
                                                                  mode = "undirected"),
                              weights = NA)
    MST_DIAM=diameter/(ncol(MST)-1) #normalised
    NTW_measures[beh_data$participant_id[sub],
                 c('MST_LEAF','MST_DIAM')]=c(MST_LEAF,MST_DIAM)
    
  }
  
  write.csv(NTW_measures,file = paste0('HCP/taskFC/',taskFC,'_NTW_measures_HCP.csv'))

}