#This script will test the effect of each topological network measure on every cognitive test score, separately this time 

library(dplyr) 

######################################################################
######standardising subjects across fMRI task and cognitive tests#####
######################################################################

#identify subjects that have data for all fMRI tasks in the 3T scan
tasklist=c('REST3T','MOTOR','GAMBLING','LG','WM','SOCIAL','RELATIONAL')
for (taskFC in tasklist)
{
 dataset_lm_resp=readRDS(paste0('datasets_model_data/HCPmodel_data_',taskFC,'.rds'))  
 #assign
 assign(paste0('subs_',taskFC), as.character(unique(dataset_lm_resp$participant_id)))
}
common_subs <- Reduce(intersect, list(subs_GAMBLING, subs_LG, subs_MOTOR, subs_RELATIONAL, subs_REST3T, subs_SOCIAL, subs_WM))
#identify subject that have data for all cognitive tests
toexclude=c()
for (i in 1:nrow(dataset_lm_resp))
{ 
  if(all(!is.na(dataset_lm_resp[i,]))==FALSE)
  {toexclude=c(toexclude,as.character(dataset_lm_resp$participant_id[i]))}
}
common_subs3T=common_subs[which(! common_subs %in% toexclude)]

#same process for subjects that had a 7T scan
tasklist=c('REST7T','MOVIE')
for (taskFC in tasklist)
{
  dataset_lm_resp=readRDS(paste0('datasets_model_data/HCPmodel_data_',taskFC,'.rds'))
  #assign
  assign(paste0('subs_',taskFC), as.character(unique(dataset_lm_resp$participant_id)))
}
common_subs <- Reduce(intersect, list(subs_REST7T, subs_MOVIE))
#identify subject that have data for all cognitive tests
toexclude=c()
for (i in 1:nrow(dataset_lm_resp))
{ 
  if(all(!is.na(dataset_lm_resp[i,]))==FALSE)
  {toexclude=c(toexclude,as.character(dataset_lm_resp$participant_id[i]))}
}
common_subs7T=common_subs[which(! common_subs %in% toexclude)]

######################################################################
####################function for assumptions check####################
######################################################################

library(car)
library(performance)
#check multicollinearity: 
#VIF values: below 10 is good, see Myers, R.H. (1990) Classical and Modern Regression Application (2nd edn). Boston: Duxbury Press.
#"tolerance" value (1/VIF): below 0.2 is a problem, cf Menard, S. (1995). Applied logistic regression analysis. Sage university paper series on quantitative applications in the social sciences, 07-106. Thousand Oaks, CA: Sage.
#Intercorrelation above 0.9 : Mayers, A. (2013). Introduction to statistics and SPSS in psychology. Edinburgh Gate, Harlow: Pearson Education Limited; Field, A. (2013). Discovering statistics using IBM SPSS statistics.).

as_check=function(dataset,model,meas){
  #RAW_FC usually is a source of multicollinearity 
  #not applicable for models with only RAW_FC, and will not check for correlation
  #specifically with RAW_FC if RAW_FC was dropped, only VIF/tolerance
  if (meas!='RAW_FC' & "RAW_FC" %in% attr(model$terms, "term.labels")) 
  {
    #check multicollinearity
    if(!all(check_collinearity(model)$VIF <= 10) | !all(1/check_collinearity(model)$VIF >= 0.2) | cor.test(dataset[[meas]],dataset[['RAW_FC']])$estimate > 0.9) 
      {
        #report indicators of multicollinearity
        cat('Multicollinearity detected\n VIF - ')
        cat(check_collinearity(model)[grep(paste0(meas,'|RAW_FC'),
                                           check_collinearity(model)$Term),'VIF'])
        cat('; Tolerance - ')
        cat(1/check_collinearity(model)[grep(paste0(meas,'|RAW_FC'),
                                             check_collinearity(model)$Term),'VIF'])
        cat(paste0('; Correlation between ',meas,' and RAW_FC - ', 
                   round(cor.test(dataset[,meas],dataset[,'RAW_FC'])$estimate,2), '\n'))
        multicollinear='yes'
     }
  } else {
    if(!all(check_collinearity(model)$VIF <= 10) | !all(1/check_collinearity(model)$VIF >= 0.2)) 
    {
      #report indicators of multicollinearity (no RAW FC)
      cat('Multicollinearity detected\n VIF - ')
      cat(check_collinearity(model)[,'VIF'])
      cat('; Tolerance - ')
      cat(1/check_collinearity(model)[,'VIF'])
      multicollinear='yes'
    }
  }
  
  #check independence of errors (Durbin-Watson statistics between 1 and 3 is good according to Mayers, A. (2013))
  DWstat=durbinWatsonTest(residuals(model))
  if (DWstat < 1 | DWstat > 3) {autocorrelation='yes'}
  
  #report violations
  violations=c()
  if (exists('multicollinear')) 
    {  
    return(c(violations, 'Assumption checks: multicollinearity found\n\n')) 
    } 
    else if (exists('autocorrelation')) 
    { 
    cat('\nAutocorrelation detected\n')
    return(c(violations, 'Assumption checks: residuals are not independent\n\n'))
    } 
  else 
    {return('\nAssumption checks: valid\n')}
}

######################################################################
###############running all models across fMRI tasks###################
######################################################################

for (taskFC in c('REST3T','REST7T','MOTOR','GAMBLING','LG','WM','SOCIAL','RELATIONAL','MOVIE'))
{
  
  cat(paste('\n######\n\nComputing models for',taskFC,'...\n\n#####\n\n'))
  
  #load network data (just to get a list of topology measures)
  NTW_measures=read.csv(paste0('HCP/taskFC/',taskFC,'_NTW_measures_HCP.csv'),row.names = 1)
  #predictors and outcome variables model data formatted in s5
  dataset_lm_resp=readRDS(paste0('datasets_model_data/HCPmodel_data_',taskFC,'.rds'))  
  
  #force all subs to have data across fMRI tasks and tests (separately for 3T/7T)
  if (taskFC=='MOVIE' | taskFC=='REST7T')
  {
    dataset_lm_resp=dataset_lm_resp[which(dataset_lm_resp$participant_id %in% common_subs7T),]
  }else{
    dataset_lm_resp=dataset_lm_resp[which(dataset_lm_resp$participant_id %in% common_subs3T),]
  }
  
  ###################################################################
  ############################Model preparation######################
  
  #standardize numerical predictors so no need to adjust beta later (better now than later: failures from parameters::standardize_parameters)
  for (u in unique(dataset_lm_resp$cog_name)) {
    dataset_lm_resp[[u]] <- as.numeric(scale(dataset_lm_resp[[u]]))}
  for (v in colnames(NTW_measures)) {
    dataset_lm_resp[[v]] <- as.numeric(scale(dataset_lm_resp[[v]]))}
  
  #function to collate statistics to report from all models across loops
  mod_tables=function(dataset_lm_resp){ 
    #prepare pvals table for correction
    p_vals=data.frame(matrix(ncol=length(unique(dataset_lm_resp$cog_name)), nrow = 1))
    colnames(p_vals)=unique(dataset_lm_resp$cog_name); row.names(p_vals)='p_vals'
    #prepare rsq/beta tables for plotting
    adjb_vals=data.frame(matrix(ncol=length(unique(dataset_lm_resp$cog_name)), nrow = 1))
    colnames(adjb_vals)=unique(dataset_lm_resp$cog_name); row.names(adjb_vals)='adjb_vals'
    rsq_vals=data.frame(matrix(ncol=length(unique(dataset_lm_resp$cog_name)), nrow = 1))
    colnames(rsq_vals)=unique(dataset_lm_resp$cog_name); row.names(rsq_vals)='rsq_vals'
    #count eligible subjects (should be consistent but still saved just in case)
    Nsubs=data.frame(matrix(ncol=length(unique(dataset_lm_resp$cog_name)), nrow = 1))
    colnames(Nsubs)=unique(dataset_lm_resp$cog_name); row.names(Nsubs)='Nsubs'
    return(rbind(p_vals, adjb_vals, rsq_vals, Nsubs))
  }
  
  #prepare global tables for plots
  val_tables_list=list()
  val_tables=mod_tables(dataset_lm_resp)
  
  ###################################################################
  ##############################Running models#######################
  
  #for each topology measure, test each cognitive score
  for (ntw_meas in colnames(NTW_measures))  
  {
    for (cog_name in unique(dataset_lm_resp$cog_name))
    {
        #narrow down dataset to selected cognitive test
        dataset_lm_resp_trial_specific=dataset_lm_resp[which(dataset_lm_resp$cog_name==cog_name),]
        
        cat(paste0('\n', ntw_meas, ' and ', cog_name,' (lm):\n'))
        
        if (ntw_meas!='RAW_FC') #as long as not RAW_FC, we additionally account for RAW_FC in the model
        {
          model=lm(paste0('performance ~ age + sex + ', ntw_meas, '+ RAW_FC'), 
                   data = dataset_lm_resp_trial_specific)
          #check assumptions
          assumption=as_check(dataset_lm_resp_trial_specific,
                              model,ntw_meas)
          if (length(grep('multicollinear',assumption))==1)
          { 
            cat(paste('Multicollinearity found, dropping RAW_FC...'))
            model=lm(paste0('performance ~ age + sex + ', ntw_meas),
                     data = dataset_lm_resp_trial_specific)
          }
        } else {
          #Modelling for RAW_FC on its own
          model=lm(paste0('performance ~ age + sex + ', ntw_meas), 
                   data = dataset_lm_resp_trial_specific)
        }
        
        #final assumption check (whether RAW_FC was dropped or not)
        cat(as_check(dataset_lm_resp_trial_specific,model,ntw_meas))
        
        #model summary
        summod=summary(model)
        #get partial rsquare
        partialrsq=suppressWarnings(unlist(rsq::rsq.partial(model)[3])[which(unlist(rsq::rsq.partial(model)[2])==ntw_meas)])
        
        print(paste0(
          cog_name,' was predicted by ', ntw_meas, 
          ' at p=', round(as.numeric(format(summod$coefficients[ntw_meas,"Pr(>|t|)"])),4),
          ', β=',  as.numeric(summod$coefficients[ntw_meas,"Estimate"]),
          ', partial R2=', round(as.numeric(partialrsq),4)
        ))
        
        #count applicable subs for this fMRI task
        Nsubs=length(unique(dataset_lm_resp_trial_specific$participant_id[which(!is.na(dataset_lm_resp_trial_specific$performance))]))
        
        #store all p values and r values across models for plotting
        val_tables[,cog_name]=
          c(summod$coefficients[ntw_meas,"Pr(>|t|)"],
            as.numeric(summod$coefficients[ntw_meas,"Estimate"]),
            as.numeric(partialrsq),
            Nsubs
            )
    }
    
    #add to global list of model tables
    val_tables_list[[ntw_meas]]=as.matrix(val_tables)
  }
  
  #Reframe all concatenated val tables in a concise table
  combined_table <- do.call(rbind, lapply(names(val_tables_list), function(meas) {
    plot_table <- as.data.frame(t(val_tables_list[[meas]]))  
    plot_table$cog_name <- row.names(plot_table)        
    plot_table$category <- meas
    plot_table$taskFC <- taskFC    
    return(plot_table)
  }))
  
  #store in broader data frame for plotting
  if(exists('general_combined_table'))
  {general_combined_table=rbind(general_combined_table,combined_table)} else
  {general_combined_table=combined_table}

}

###############################################################
###############################################################
##############################PLOTTING#########################

library(ggplot2)
combined_table=general_combined_table 

#factorise type of cognitive test
combined_table$cog_name <- factor(combined_table$cog_name, levels = unique(combined_table$cog_name))
#remove potential negative rsquares (possible modelling errors due to it being extremely low)
combined_table[combined_table[, "rsq_vals"] < 0, "rsq_vals"] <- 0

#assign name of topology measures' categories for colouring in plot
combined_table$ntw_type <- NA
for (ntnum in 1:length(combined_table$category)) {
  nt=combined_table$category[ntnum]
  if (nt=='CLUSTERING0.1'|nt=='CLUSTERING0.2'|nt=='CLUSTERING0.3'|nt=='GLOB_EFF0.1'|nt=='GLOB_EFF0.2'|nt=='GLOB_EFF0.3')
    {combined_table$ntw_type[ntnum]='Graph measures'}
  if (nt=='MST_DIAM'|nt=='MST_LEAF')
    {combined_table$ntw_type[ntnum]='Minimum Spanning Tree'}
  if (nt=='PH_BD'|nt=='PH_BS'|nt=='PH_CS')
    {combined_table$ntw_type[ntnum]='Persistent Homology'}
  if (nt=='RAW_FC')
    {combined_table$ntw_type[ntnum]='Raw functional connectivity'}
}

#rename network measures for clarity in the plot
category_names <- c(
  "GLOB_EFF0.1" = "Global efficiency (t.10%)",
  "CLUSTERING0.1" = "Clustering coefficient (t.10%)",
  "GLOB_EFF0.2" = "Global efficiency (t.20%)",
  "CLUSTERING0.2" = "Clustering coefficient (t.20%)",
  "GLOB_EFF0.3" = "Global efficiency (t.30%)",
  "CLUSTERING0.3" = "Clustering coefficient (t.30%)",
  "PH_BS" = "Backbone Strength",
  "PH_CS" = "Cycle Strength",
  "PH_BD" = "Backbone Dispersion",
  "MST_LEAF" = "Leaf fraction",
  "MST_DIAM" = "Diameter",
  "RAW_FC" = "Mean connectivity"
)
# Replace abbreviations with full names in the `category` column
combined_table$category_abbr <- combined_table$category
combined_table$category <- category_names[as.character(combined_table$category_abbr)]

#rename fMRI tasks for clarity in the plot
combined_table$taskFC[which(combined_table$taskFC=='LG')]='LANGUAGE'
combined_table$taskFC[which(combined_table$taskFC=='WM')]='WORKING\nMEMORY'
combined_table$taskFC[which(combined_table$taskFC=='SOCIAL')]='SOCIAL\n COGNITION'
combined_table$taskFC[which(combined_table$taskFC=='RELATIONAL')]='RELATIONAL\nPROCESSING'
combined_table$taskFC[which(combined_table$taskFC=='REST3T')]='REST (3T)'
combined_table$taskFC[which(combined_table$taskFC=='REST7T')]='REST (7T)'
combined_table$taskFC <- factor(combined_table$taskFC,
levels = c("GAMBLING", "LANGUAGE", "MOTOR", "RELATIONAL\nPROCESSING", "SOCIAL\n COGNITION", "WORKING\nMEMORY", "REST (3T)", "REST (7T)", "MOVIE"))

#categorise per topology measure and fMRI task
combined_table$task_category <- paste0(combined_table$category, "_", 
                                       combined_table$taskFC)
combined_table <- combined_table %>%
  mutate(task_category = factor(task_category,
  levels = unique(task_category[order(taskFC, ntw_type, category)])))
#also create grouping and ordering for 3T and 7T tasks distinctly
combined_table$task_group <- ifelse(combined_table$taskFC %in% c("REST (7T)", "MOVIE"), "Rest/Movie","Other Tasks")
#order accordingly
combined_table <- combined_table %>%
  mutate(task_category = factor(task_category,
  levels = unique(task_category[order(task_group, taskFC, ntw_type, category)])))

# Get x-axis positions for each fMRI task block
taskFC_breaks <- combined_table %>%
  group_by(taskFC) %>%
  summarise(x_break = max(as.numeric(task_category)) + 0.5) %>%
  arrange(x_break)
# Find the boundary between 3T tasks and the 7T tasks (REST 7T first to the left)
rest7t_left_edge <- combined_table %>%
  filter(taskFC == "REST (7T)") %>%
  summarise(x_pos = min(as.numeric(task_category))) %>%
  pull(x_pos)

#get max rsquare of each dimension for later computation, also ymin for plot
max_rsq_per_facet <- combined_table %>%
  group_by(cog_name, taskFC) %>%
  summarise(max_rsq = max(rsq_vals, na.rm = TRUE), .groups = "drop")
min_rsq_per_facet <- combined_table %>%
  group_by(cog_name, taskFC) %>%
  summarise(min_rsq = min(rsq_vals, na.rm = TRUE), .groups = "drop")
combined_table <- combined_table %>%
  left_join(max_rsq_per_facet, by = c("cog_name", "taskFC")) %>%
  left_join(min_rsq_per_facet, by = c("cog_name", "taskFC"))

#this adds the offset for relative position of beta value
combined_table <- combined_table %>%
  group_by(cog_name) %>%
  mutate(y_offset = (max(rsq_vals) - min(rsq_vals)) * 0.1)  
combined_table <- combined_table %>% mutate(label_hjust = 0.5) 

#add FDR corrected p-values
combined_table$p_vals_fdr=p.adjust(combined_table$p_vals,method = 'fdr')

# Plot with corrected grouping
p3 <- ggplot(combined_table, aes(x = task_category, y = rsq_vals, group = ntw_type)) +
  geom_ribbon(aes(x=as.numeric(task_category), ymin = min_rsq, ymax = rsq_vals, 
                  group = interaction(taskFC, ntw_type), fill = ntw_type), 
              alpha = 0.2) +  # Ensure lines match header colors
  geom_line(aes(group = interaction(taskFC, ntw_type), color = ntw_type), linewidth = 0.7) + #add thick line matching category colors
  # Non-significant points correctly colored
  geom_point(aes(fill = ntw_type, shape = ntw_type, color = ntw_type), size = 3, stroke = 0.5) +
  geom_point(data = subset(combined_table, p_vals_fdr < 0.05), 
             aes(fill = ntw_type), shape = 8, size = 2, stroke = 1, 
             color = "black", show.legend = FALSE) + #Significant asterisks, FDR
  geom_point(data = subset(combined_table, p_vals < 0.05), 
             aes(fill = ntw_type), shape = 8, size = 0.5, stroke = 1, 
             color = "black", show.legend = FALSE) + #not FDR corrected: smaller
  geom_text(aes(label = paste0('β=', round(adjb_vals, 2)),
                y = rsq_vals + y_offset,
                hjust = label_hjust),
            size = 2, color = "#5E5E5E",
            angle=90, vjust=0.5) +   #beta values
  labs(
    x = "Network Measure",
    y = "Partial R-Squared",
    fill = "Category",
    shape = "Category",
    color = "Category"
  ) +
  facet_wrap(~cog_name, ncol = 3, scales = "free_y",
             labeller = label_wrap_gen(width = Inf)) + #Y scale relative to each panel
  theme_minimal() +
  theme(plot.margin = margin(10, 10, 10, 20))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.text = element_text(size = 10, face = "bold", 
                                  lineheight = 0.2)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_fill_manual(values = c(
    "Graph measures" = "#3468A4",
    "Persistent Homology" = "#95435C",
    "Minimum Spanning Tree" = "#C0915C",
    "Raw functional connectivity" = "#9933FF"
  )) +
  scale_shape_manual(values = c(
    "Graph measures" = 16,
    "Persistent Homology" = 15,
    "Minimum Spanning Tree" = 18,
    "Raw functional connectivity" = 17
  )) +
  scale_color_manual(values = c(
    "Graph measures" = "#3468A4",
    "Persistent Homology" = "#95435C",
    "Minimum Spanning Tree" = "#C0915C",
    "Raw functional connectivity" = "#9933FF"
  )) +
  theme(
    legend.position = "bottom", 
    legend.justification = "centre",
    legend.margin = margin(-10, 0, 0, 0)) + 
  guides(fill = guide_legend(override.aes = list(stroke = 0)))  # Keep legend clean, no halos

#add dashed lines between taskFC
p3 <- p3 +
  geom_vline(data = taskFC_breaks,
             aes(xintercept = x_break),
             linetype = "dashed", color = "grey50", linewidth = 0.5)
# Add thick line between REST7T and MOVIE
p3 <- p3 +
  geom_vline(xintercept = rest7t_left_edge - 0.5, color = "grey50", linewidth = 1)

#print taskFC name between the dashed lines
facet_max_y <- combined_table %>%
  group_by(cog_name) %>%
  summarise(y_pos = max(rsq_vals, na.rm = TRUE), .groups = "drop")
taskFC_labels <- combined_table %>%
  group_by(cog_name, taskFC) %>%
  summarise(x_pos = mean(as.numeric(task_category)), .groups = "drop") %>%
  left_join(facet_max_y, by = "cog_name")
p3 <- p3 +
  geom_text(data = taskFC_labels,
            aes(x = x_pos, y = y_pos, label = taskFC),
            inherit.aes = FALSE,
            color = "grey30", alpha = 0.15,
            size = 6, fontface = "bold",
            family = "Impact")  

#don't need fMRI task name in the plot's x labels as named inside the facet
label_lookup <- combined_table %>% select(task_category, category) %>% distinct()
p3 <- p3 + scale_x_discrete(labels = setNames(label_lookup$category, label_lookup$task_category))
print(p3)

#split in 3 for clarity
source('#jobs/s7b_splitter.R')

#whole unsplitted image
ggsave(paste0("figures/taskFC/HCP_task-specific_perf_alltasks.png"), 
       plot = p3, units = 'px', width=15000, height=6000, limitsize=F)

#save .rds for further analyses
saveRDS(combined_table,file=paste0('datasets_stat/taskFC/stat_output_HCP_alltasks.rds'))
