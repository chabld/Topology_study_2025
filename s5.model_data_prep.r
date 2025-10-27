#This script formats and collates topology measures (from s4) and behavioural data appropriately before statistics are run on them
#requires all samples to have the [taskFC]_NTW_measures_HCP.csv as per s4 and the behavioural data HCP_behdata.csv


for (taskFC in c('REST3T','REST7T','MOTOR','GAMBLING','LG','WM','SOCIAL','RELATIONAL','MOVIE'))
{
  #get topology measures precomputed in s4
  NTW_measures=read.csv(paste0('HCP/taskFC/',taskFC,'_NTW_measures_HCP.csv'),row.names = 1)
  
  #get matching participant data
  beh_data <- read.csv(paste0('HCP/HCP_behdata.csv'),row.names = 1)
  colnames(beh_data)[1]='participant_id'
  beh_data=beh_data[which(beh_data$participant_id %in% row.names(NTW_measures)),] 
  
  #select cognitive measures
  meas_of_interest=c('PicSeq_Unadj', 'CardSort_Unadj',
                     'Flanker_Unadj', 'PMAT24_A_CR',
                     'ReadEng_Unadj','PicVocab_Unadj',
                     'ProcSpeed_Unadj', 'DDisc_AUC_200',
                     'VSPLOT_TC','SCPT_TP', 'IWRD_TOT',
                     'ListSort_Unadj','CogFluidComp_Unadj',
                     'CogCrystalComp_Unadj', 'ER40_CR')
  
  #collate model variables
  dataset_lm=cbind(participant_id=as.factor(beh_data$participant_id),
                   age=as.factor(beh_data[,'Age']), #HCP's age are age-bins
                   sex=as.factor(beh_data[,'Gender']),
                   beh_data[,meas_of_interest],
                   NTW_measures)
  
  #rename cognitive tests for clarity
  meas_of_interest_renamed=c('Episodic Memory (Picture Sequence Memory)', 'Cognitive Flexibility (Dimensional Change Card Sort)', 'Inhibition (Flanker Task)', 'Fluid Intelligence (Penn Progressive Matrices)', 'Oral Reading Recognition', 'Vocabulary Comprehension  (Picture Vocabulary)','Processing Speed  (Pattern Completion Processing Speed)', 'Self-regulation/Impulsivity (Delay Discounting 200k)', 'Spatial Orientation (Variable Short Penn Line Orientation Test)', 'Sustained Attention (Short Penn Continuous Performance Test)', 'Verbal Episodic Memory (Penn Word Memory Test)', 'Working Memory', 'Cognition Fluid Composite', 'Cognition Crystallized Composite', 'Emotion Recognition (Penn Emotion Recognition Test)')
  colnames(dataset_lm)[colnames(dataset_lm) %in% meas_of_interest] = meas_of_interest_renamed
  
  #all cognitive measures stacked in columns
  cog_stacked=stack(dataset_lm[,meas_of_interest_renamed])
  cog_stacked[,2]=as.factor(cog_stacked[,2])
  colnames(cog_stacked)=c('performance', 'cog_name')
  dataset_lm_resp=cbind(dataset_lm, cog_stacked)
  
  saveRDS(dataset_lm_resp, paste0('datasets_model_data/HCPmodel_data_',taskFC,'.rds'))

}
