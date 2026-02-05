#This script will test effect of each topological network measure on scores of all cognitive measures merged together

library(dplyr)
library(nlme)

######################################################################
#####standardising subjects across fMRI tasks and cognitive tests#####
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
                   round(cor.test(dataset[[meas]],dataset[['RAW_FC']])$estimate,2), '\n'))
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
      stop()
    }
  }
  #check independence of errors (Durbin-Watson statistics between 1 and 3 is good according to Mayers, A. (2013))
  DWstat=durbinWatsonTest(as.numeric(residuals(model, type = "normalized"))) #adapted to lme
  if (DWstat < 1 | DWstat > 3) {autocorrelation=TRUE}
  
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
################running all models across fMRI tasks##################
######################################################################

#prepare global table for plot
val_tables_list=list()

for (taskFC in c('REST3T','REST7T','MOTOR','GAMBLING','LG','WM','SOCIAL','RELATIONAL','MOVIE'))
{
  cat(paste('\n######\n\nComputing models for',taskFC,'fMRI...\n\n#####\n\n'))
  
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
  
  ##############################################################
  #####################Model preparation########################
  
  #turn each performance to a Z score for interpretability across measures
  dataset_lm_resp$performance <- dataset_lm_resp %>%
    group_by(cog_name) %>%
    mutate(z = as.numeric(scale(performance))) %>%
    pull(z)
  #standardize numerical predictors so no need to adjust beta later (better now than later: failures from parameters::standardize_parameters)
  for (v in colnames(NTW_measures)) {
    dataset_lm_resp[[v]] <- as.numeric(scale(dataset_lm_resp[[v]]))}
  
  #function to collate statistics to report from all models across loops
  mod_tables=function(NTWmeas){ 
    #prepare pvals table for correction
    p_vals=data.frame(matrix(ncol=length(NTWmeas), nrow = 1))
    colnames(p_vals)=colnames(NTWmeas); row.names(p_vals)='p_vals'
    #prepare rsq/beta tables for plotting
    adjb_vals=data.frame(matrix(ncol=length(NTWmeas), nrow = 1))
    colnames(adjb_vals)=colnames(NTWmeas); row.names(adjb_vals)='adjb_vals'
    rsq_vals=data.frame(matrix(ncol=length(NTWmeas), nrow = 1))
    colnames(rsq_vals)=colnames(NTWmeas); row.names(rsq_vals)='rsq_vals'
    #prepare Nsubs tables for plotting
    Nsubs=data.frame(matrix(ncol=length(NTWmeas), nrow = 1))
    colnames(Nsubs)=colnames(NTWmeas); row.names(Nsubs)='Nsubs'
    return(rbind(p_vals, adjb_vals, rsq_vals, Nsubs))
  }

  #prepare tables
  val_tables=mod_tables(NTW_measures)
  mod_list_rest=list()
  
  ###################################################################
  ##############################Running models#######################
  
  #safety removal if NA as nlme's lme() does not handle them
  dataset_lm_resp=dataset_lm_resp[which(!is.na(dataset_lm_resp$performance)),]
  
  #homogeneity of variance was violated for cog_name, so added varIdent in lme()
  #levenetest=leveneTest(as.formula(paste0('performance ~ cog_name')), data = dataset_lm_resp)
  #if(levenetest$`Pr(>F)`[1] < .05) {cat('Heteroscedasticity: cog_name\'s levels do not have homogeneous variance. Will account for that in the models.\n\n')}
  
  for (ntw_meas in colnames(NTW_measures))
  {
    if (ntw_meas!='RAW_FC') #as long as not RAW_FC, we additionally account for RAW_FC in the model
    {
        tryCatch({ #had to catch error here as RAW_FC's collinearity made some models crash before assumption checks could be done
          model <- nlme::lme(
            as.formula(paste0('performance ~ age + sex + cog_name +', 
                              ntw_meas,' + RAW_FC')), 
            random = ~1 | participant_id,
            #corCompSymm fixes autocorrelation of residuals 
            correlation = corCompSymm(form = ~1 | participant_id),
            control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
            #homoscedasticity not assumed — residuals may have different variance between cognitive tests
            weights = varIdent(form = ~1 | cog_name),
            data = dataset_lm_resp)
        }, error = function(e) 
        {
          message(paste("Not including RAW_FC due to convergence error for", ntw_meas, "— fitting reduced model.")); 
          assumption <<- 'fail'
        })
        
        #check assumptions but only drop RAW_FC if failure/multicollinearity happened
        if (exists('assumption')) 
        {if(assumption!='fail') {assumption=as_check(dataset_lm_resp,model,ntw_meas)}
        }else{assumption=as_check(dataset_lm_resp,model,ntw_meas)}
        #if multicollinearity or model failed due to RAW_FC, drop RAW_FC
        if (length(grep('multicollinear',assumption))==1|assumption=='fail')
        { 
          cat(paste('Multicollinearity found, dropping RAW_FC...'))
          tryCatch({ #had to catch error as may STILL struggle if minuscule effect
            model <- nlme::lme(
              as.formula(paste0('performance ~ age + sex + cog_name + ',
                                ntw_meas)), 
              random = ~1 | participant_id,
              correlation = corCompSymm(form = ~1 | participant_id),
              control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
              weights = varIdent(form = ~1 | cog_name),
              data = dataset_lm_resp)
          }, error = function(e) 
          {
            message(paste("Not including RAW FC but still failed. Trying to omit corCompSymm adjustment...")); 
            assumption <<- 'fail'
          })
          
          if (assumption=='fail'){
            model <- nlme::lme(
              as.formula(paste0('performance ~ age + sex + cog_name + ',
                                ntw_meas)), 
              random = ~1 | participant_id,
              control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
              weights = varIdent(form = ~1 | cog_name),
              data = dataset_lm_resp)
          }
        }
    } else {
      #Modelling for RAW_FC on its own
      model <- tryCatch({ #similar safety in case of crash
        nlme::lme(
            as.formula(paste0('performance ~ age + sex + cog_name +', ntw_meas)), 
            random = ~1 | participant_id,
            correlation = corCompSymm(form = ~1 | participant_id),
            control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
            weights = varIdent(form = ~1 | cog_name),
            data = dataset_lm_resp)
        }, error = function(e) {
          #In case we want to tweak parameters if it fails. But no crash was observed
          #with our current data so the stayed the same.
          nlme::lme(
            as.formula(paste0('performance ~ age + sex + cog_name +', ntw_meas)),             random = ~1 | participant_id,
            correlation = corCompSymm(form = ~1 | participant_id),
            control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
            weights = varIdent(form = ~1 | cog_name),
            data = dataset_lm_resp)
        })
   }
  
  #final assumption check (whether RAW_FC was dropped or not)
  cat(as_check(dataset_lm_resp,model,ntw_meas))
  
  #model summary
  summod=summary(model)
  partialrsq=r2glmm::r2beta(model)
  partialrsq=partialrsq[which(partialrsq$Effect==ntw_meas),'Rsq']
  
  print(paste0(
    'Performance was predicted by ', ntw_meas, 
    ' at p=', round(as.numeric(format(summod$tTable[ntw_meas,"p-value"])),4),
    ', β=', as.numeric(summod$tTable[ntw_meas,"Value"]),
    ', partial R2=', round(as.numeric(partialrsq),4)
  ))
  
  #count applicable subs for this fMRI task
  Nsubs=length(unique(dataset_lm_resp$participant_id[which(!is.na(dataset_lm_resp$performance))]))
  
  #store all p values and r values across models for plotting
  val_tables[,ntw_meas]=
    c(summod$tTable[ntw_meas,"p-value"],
      as.numeric(summod$tTable[ntw_meas,"Value"]),
      as.numeric(partialrsq),
      Nsubs
      )
  
  #store all models
  mod_list_rest[[ntw_meas]]=model
  }
  
  #remove columns not tested (if applicable) as they end up NA
  val_tables=val_tables[,c(names(mod_list_rest))]
  #add to global list of model tables
  val_tables_list[[taskFC]]=as.matrix(val_tables) 
}

#save .rds for potential replotting
saveRDS(val_tables_list,file=paste0('datasets_stat/HCP_all_cog_val_tables_list.rds'))

###############################################################
###############################################################
##############################PLOTTING#########################

library(ggplot2)
#Reframe all concatenated val tables in a concise table
combined_table <- do.call(rbind, lapply(names(val_tables_list), function(meas) {
  plot_table <- as.data.frame(t(val_tables_list[[meas]]))  
  plot_table$category <- row.names(plot_table)        
  plot_table$taskFC <- meas               
  return(plot_table)
}))
#factorise topology measure
combined_table$category <- factor(combined_table$category, levels = unique(combined_table$category))
combined_table$category_abbr <- combined_table$category

#assign name of topology measures' categories for colouring in plot
combined_table$ntw_type <- NA   
for (ntnum in seq_along(combined_table$category_abbr)) {
  nt <- combined_table$category_abbr[ntnum]
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
combined_table$category <- category_names[as.character(combined_table$category_abbr)]

#Order category by ntw_type so they are grouped together
library(dplyr)
ordered_levels <- combined_table %>% distinct(category, ntw_type) %>%
  arrange(ntw_type, category) %>% pull(category)
combined_table$category <- factor(combined_table$category, levels = ordered_levels)

#rename fMRI tasks for clarity in the plot
combined_table$taskFC[which(combined_table$taskFC=='LG')]='LANGUAGE'
combined_table$taskFC[which(combined_table$taskFC=='WM')]='WORKING MEMORY'
combined_table$taskFC[which(combined_table$taskFC=='SOCIAL')]='SOCIAL COGNITION'
combined_table$taskFC[which(combined_table$taskFC=='RELATIONAL')]='RELATIONAL PROCESSING'
combined_table$taskFC[which(combined_table$taskFC=='REST3T')]='RESTING-STATE (3T)'
combined_table$taskFC[which(combined_table$taskFC=='REST7T')]='RESTING-STATE (7T)'

#get max rsquare of each dimension so plot can account for that in the ribbons
max_rsq_per_facet <- combined_table %>%
  group_by(taskFC) %>%
  summarise(max_rsq = max(rsq_vals, na.rm = TRUE))
combined_table <- left_join(combined_table, max_rsq_per_facet, by = "taskFC")
min_rsq_per_facet <- combined_table %>%
  group_by(taskFC) %>%
  summarise(min_rsq = min(rsq_vals, na.rm = TRUE))
combined_table <- left_join(combined_table, min_rsq_per_facet, by = "taskFC")

#this adds the offset for relative position of beta value
global_range <- max(combined_table$rsq_vals, na.rm = TRUE) -
  min(combined_table$rsq_vals, na.rm = TRUE)
combined_table <- combined_table %>%
  mutate(y_offset = global_range * 0.1)
#Tracking half of Y axis to know if beta should be plotted above or below main value
global_mid_y <- (max(combined_table$rsq_vals, na.rm = TRUE) +
                   min(combined_table$rsq_vals, na.rm = TRUE)) / 2
combined_table <- combined_table %>%
  mutate(label_hjust = ifelse(rsq_vals > global_mid_y, 1.6, -0.1))   #1 = above half, 0.1 = below

#add FDR corrected p-values
combined_table$p_vals_fdr=p.adjust(combined_table$p_vals,method = 'fdr')

# Plot with corrected grouping
p3 <- ggplot(combined_table, aes(x = category, y = rsq_vals, group = ntw_type)) +
  geom_ribbon(aes(ymin = min_rsq, ymax = rsq_vals, group = ntw_type, fill = ntw_type),
              alpha = 0.2) + #groups topology measure categories in ribbons
  geom_line(aes(color = ntw_type), linewidth = 0.7) +  #add thick line matching category colors
  geom_point(aes(fill = ntw_type, shape = ntw_type, color = ntw_type), size = 3, stroke = 0.5) +  #ensures all points also match category color
  geom_point(data = subset(combined_table, p_vals_fdr < 0.05), 
             aes(fill = ntw_type), shape = 8, size = 3, stroke = 1, 
             color = "black", show.legend = FALSE) + #Significant asterisks, FDR
  geom_point(data = subset(combined_table, p_vals < 0.05), 
             aes(fill = ntw_type), shape = 8, size = 0.5, stroke = 1, 
             color = "black", show.legend = FALSE) + #not FDR corrected: smaller
  geom_text(aes(label = paste0('β=', round(adjb_vals, 2)),
              y = rsq_vals + y_offset,
              hjust = label_hjust),  
            size = 3, color = "#5E5E5E",
            angle=90, vjust=0.5) + #beta values
  labs(
    x = "Network topology measure",
    y = "Semi-partial R-Squared",
    fill = "Category",
    shape = "Category",
    color = "Category"
  ) +
  facet_wrap(~taskFC, ncol = 2, scales = "fixed",
             labeller = label_wrap_gen(width = Inf)) + #Y scale relative to each panel
  theme_minimal() +
  theme(plot.margin = margin(10, 10, 10, 20))+
  theme(axis.title.y = element_text(size = 14)) + 
  theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1)) +
  theme(strip.text = element_text(size = 12, face = "bold", 
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
    legend.margin = margin(-10, 45, 0, 0)) + 
  guides(fill = guide_legend(override.aes = list(stroke = 0)))  #Keeps legend clean, no halos

#Add number of subjects inside panels for specificity
global_max_y <- max(combined_table$rsq_vals, na.rm = TRUE) + 0.002
facet_labels <- combined_table %>%
  group_by(taskFC) %>%
  summarise(
    Nsubs = unique(Nsubs),
    .groups = "drop"
  ) %>%
  mutate(label_y = global_max_y)  
left_category <- levels(factor(combined_table$category))[2]
# Adds label without expanding y-axis
p3 <- p3 +
  geom_text(data = facet_labels,
            aes(x = left_category, y = label_y, label = paste0("N=", Nsubs)),
            inherit.aes = FALSE,
            size = 4, color = "grey", fontface = "italic")
print(p3)

ggsave(paste0("figures/taskFC/HCP_overall_perf_alltasks.png"), plot = p3, width = 2200, height = 2700, units = "px", dpi = 300)
