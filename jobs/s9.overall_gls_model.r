#This script tests the effect of network measures on their respective partil R-squared values obtained across all models of the study (following s7). 
#Tested separately for 3T and 7T scanners' models

library(nlme)
library(emmeans)

####################################################################
#######################dataset preparation##########################
####################################################################

for (scan in c("3T", "7T"))
{

  dataset_stat=readRDS('datasets_stat/taskFC/stat_output_HCP_alltasks.rds')
  
  if(scan=="3T")
  {
    dataset_stat=dataset_stat[which(dataset_stat$taskFC != 'REST (7T)' & dataset_stat$taskFC != 'MOVIE'),]
  } else if (scan=="7T")
  {  
    dataset_stat=dataset_stat[which(dataset_stat$taskFC == 'REST (7T)' | dataset_stat$taskFC == 'MOVIE'),]
  }
  
  #rename fMRI tasks for clarity for later plot
  levels(dataset_stat$taskFC)[levels(dataset_stat$taskFC) == "LANGUAGE"] <- "LG"
  levels(dataset_stat$taskFC)[levels(dataset_stat$taskFC) == "WORKING\nMEMORY"] <- "WM"
  levels(dataset_stat$taskFC)[levels(dataset_stat$taskFC) == "SOCIAL\n COGNITION"] <- "SOCIAL"
  levels(dataset_stat$taskFC)[levels(dataset_stat$taskFC) == "RELATIONAL\nPROCESSING"] <- "RELATIONAL"
  
  #factorise
  dataset_stat$cog_name=as.factor(dataset_stat$cog_name)
  dataset_stat$taskFC=as.factor(dataset_stat$taskFC)
  dataset_stat$taskFC <- droplevels(dataset_stat$taskFC) #drop as 3T/7T mismamtch
  dataset_stat$category=as.factor(dataset_stat$category_abbr)
  
  #assumption check function
  library(car)
  as_check=function(model){
    flags=c()
    #"type" argument added because of the interaction term
    #using adjusted gvif as more reliable when many dummy variables (categorical levels), and square root of the criterion
    #https://www.bookdown.org/rwnahhas/RMPH/mlr-collinearity.html
    if(!all(vif(model, type='predictor')[,3] <= sqrt(10))) {
      cat('Multicollinearity detected\n VIF - ')
      cat(vif(model)[,1])
      cat(';\n Tolerance - ')
      cat(1/vif(model)[,1])
      cat(';\n')
      flags=c(flags,'multicollinearity')
    }
    #check independence of errors (Durbin-Watson statistics between 1 and 3 is good according to Mayers, A. (2013))
    DWstat=durbinWatsonTest(as.numeric(residuals(model)))
    if (DWstat < 1 | DWstat > 3)  {flags=c(flags,'autocorrelated residuals')}
    
    return(flags)
  }
  
  ###############################################################
  #############################MODEL#############################
  
  #The baseline is defined as the grand mean as choosing a level as baseline is arbitrary for types of cognitive tests and of topology measures
  #This sets global defaults for how R encodes factor variables when fitting models. Specifically: contr.sum" applies sum (effect) coding to unordered factors 
  contrasts(dataset_stat$cog_name) <- contr.sum
  contrasts(dataset_stat$category) <- contr.sum
  
  #for taskFC, REST is defined as the baseline for the fMRI taks category
  if (scan=="7T")
  {
    dataset_stat$taskFC <- relevel(factor(dataset_stat$taskFC), ref = "REST (7T)")
    contrasts(dataset_stat$taskFC) <- contr.treatment(levels(dataset_stat$taskFC),
                                base = which(levels(dataset_stat$taskFC) == "REST (7T)"))
  } else 
  {     
    dataset_stat$taskFC <- relevel(factor(dataset_stat$taskFC), ref = "REST (3T)")
    contrasts(dataset_stat$taskFC) <- contr.treatment(levels(dataset_stat$taskFC),
                                base = which(levels(dataset_stat$taskFC) == "REST (3T)")) 
  }
  
  #log transformation fixes most heteroskedasticity 
  leveneTest(log(rsq_vals) ~ category, data = dataset_stat)$`Pr(>F)`[1] < .05 
  leveneTest(log(rsq_vals) ~ cog_name, data = dataset_stat)$`Pr(>F)`[1] < .05 
  leveneTest(log(rsq_vals) ~ taskFC, data = dataset_stat)$`Pr(>F)`[1] < .05 
  leveneTest(log(rsq_vals) ~ category * taskFC, data = dataset_stat)$`Pr(>F)`[1] < .05
  
  #Generalized Least Squares (GLS)

  #Because gls() refers to the source dataset_stat variable even if we assign() it to a new variable name, meta_model_3T since will refer to the replaced dataset_stat from the last loop instead of the initial object. To avoid this, we create explicitly a variable instead of assigning:
  if(scan=="3T")
  {
    dataset_stat_3T <- dataset_stat
    
    #weights controls for variance of taskFC (as it still creates heteroscedasticity after log transform)
    meta_model <- nlme::gls(
      log(rsq_vals) ~ cog_name + category * taskFC,
      data = dataset_stat_3T,
      weights = varIdent(form = ~1 | taskFC)  # Allows different variance per group
    )
    as_check(meta_model)
  } else if (scan=="7T") 
  {
    dataset_stat_7T <- dataset_stat
    
    meta_model <- nlme::gls(
      log(rsq_vals) ~ cog_name + category * taskFC,
      data = dataset_stat_7T,
      weights = varIdent(form = ~1 | taskFC)  # Allows different variance per group
    )
    as_check(meta_model) 
 }
  
  #get summary of model contrasts (estimates, p values)
  #DevContrStats script extracts the categorical levels omitted from the default summary and returns more readable results (the last last levels are hidden in R due to contr.sum and simply implied to build a contrast matrix summing to 0, see for more context:
  #https://stackoverflow.com/questions/72820236/comparing-all-factor-levels-to-the-grand-mean-can-i-tweak-contrasts-in-linear-m)
  source("#jobs/DevContrStats_gls.r")
  
  #cog task and network measure main effects
  coefficients=DevContrStats_gls(dataset_stat, meta_model, 'cog_name')
  coefficients=rbind(coefficients,
               DevContrStats_gls(dataset_stat, meta_model, 'category'))
  
  #script not needed for taskFC variable as there are only 2 levels (contr.treatment)
  #taskFC levels first  
  summod=as.data.frame(summary(meta_model)$tTable)
  coefficients=rbind(coefficients,
                     summod[grep('^taskFC',row.names(summod)),])
  
  #interactions terms 
  coefficients=rbind(coefficients,
                     DevContrStats_gls(dataset_stat, meta_model, 'category','taskFC'))
  #taskFC REST is dropped
  coefficients=coefficients[!grepl(':REST',row.names(coefficients)),]

  #record which scan each contrast is based on
  coefficients$scan=scan
  #save summaries separately per scanners
  assign(paste0('coefficients_',scan), coefficients)
  assign(paste0('meta_model_',scan), meta_model)
}

#merged again for the sake of FDR correction
coefficients_all=rbind(coefficients_3T, coefficients_7T)
#no need to correct across intercepts so ignore their coefficient
coefficients_all=coefficients_all[-grep('Intercept',row.names(coefficients_all)),]
#FDR correction
coefficients_all$`p-value_fdr`=p.adjust(coefficients_all[,"p-value"], method='fdr')
#filter out insignificant contrasts
sig_coefficients=coefficients_all[which(coefficients_all$`p-value_fdr`<.05),]
print(sig_coefficients)

###############################################################
#####################posthoc contrast##########################
###############################################################

for (scan in c("3T", "7T"))
{
  if (scan=='3T'){dataset=dataset_stat_3T; model=meta_model_3T}
  if (scan=='7T'){dataset=dataset_stat_7T; model=meta_model_7T}
    
  #test pairwise contrasts with emmeans 
  #warning about contrast levels dropped is normal, due to the contr.sum parameter
  dataset <- droplevels(dataset)
  #mode = df.error to keep same df throughout instead of computing a new one per contrast (leads to near 0 dfs which results in absurd SEs in this data), more stable for gls()
  EMM <- emmeans::emmeans(model, ~ category * taskFC, data = dataset,
                 mode='df.error')
  pairwise_comp <- emmeans::pairs(EMM, infer = TRUE, adjust = "none")
  
  #don't count pairwise comparisons within own topology measure category:
  library(stringr)
  contrasts=as.data.frame(pairwise_comp)
  duplicates=which(str_count(contrasts$contrast,'Minimum|Leaf|Diameter|LEAF|DIAM')==2 | 
                   str_count(contrasts$contrast,'Persistent|Backbone|Cycle|BS|BD|CS')==2 |
                   str_count(contrasts$contrast,'Clustering|Global|GLOB|CLUST')==2)
  if (length(duplicates)!=0){contrasts=contrasts[-duplicates,]}
  
  #specify scanner it came from
  contrasts$scan=scan
  assign(paste0('contrasts_',scan),contrasts)
}

#merge for FDR correction
contrasts=rbind(contrasts_3T, contrasts_7T)
#significant pre and post fdr correction:
contrasts$p.value_fdr=p.adjust(contrasts$p.value, method='fdr')
contrasts_uncor=contrasts
#apply correction
contrasts=contrasts[which(contrasts$p.value_fdr<.05),]
#then assign fdr corrected p back to each respective set of contrasts
contrasts_3T=contrasts[which(contrasts$scan=='3T'),]
contrasts_7T=contrasts[which(contrasts$scan=='7T'),]

###############################################################
#############################Plots#############################

library(ggplot2)
library(dplyr)
library(forcats)
library(stringr)

for (scan in c("3T", "7T"))
{

  contrasts=get(paste0("contrasts_", scan))
  #identify first and second variable in each contrast
  df_contrasts <- as.data.frame(contrasts) %>%
    mutate(original_var1 = str_trim(str_extract(contrast, "^[^-]+")),
           original_var2 = str_trim(str_extract(contrast, "(?<= - ).*")) ) %>% 
  #order in terms of increases
  mutate( flip = estimate < 0, 
          contrast = if_else(flip, paste(original_var2, "-", original_var1), contrast), 
          estimate = if_else(flip, -estimate, estimate), 
          t.ratio = if_else(flip, -t.ratio, t.ratio) )
  
  #add asterisks for significance
  df_contrasts <- df_contrasts %>%
    mutate(sig = case_when(
      p.value_fdr < 0.001 ~ "***",
      p.value_fdr < 0.01  ~ "**",
      p.value_fdr < 0.05  ~ "*",
      TRUE                ~ "" #catch-all if nothing above is TRUE
    )) 
  
  #clarify direction of difference (- replaced with <>)
  df_contrasts <- df_contrasts %>% 
    mutate(
      contrast_label = str_replace(contrast, "(.*) - (.*)", "\\1 <> \\2"), 
      contrast_label = if_else(estimate > 0, str_replace(contrast_label, "<>", " > "), str_replace(contrast_label, "<>", " < "))) %>% 
    mutate(contrast_label = fct_inorder(contrast_label)) 
  df_contrasts <- df_contrasts %>% mutate(label_leader = str_trim(str_extract(contrast_label, "^[^>]+")) # gets left-hand side of label 
          )
  df_contrasts <- df_contrasts %>% arrange(label_leader)
  
  #Colour contrast depending on strongest measure in each comparison contrast
  #Accounts for different possible labels
  color_map <- c(
    "Graph measures" = "#3468A4",
    "Clustering coefficient" = "#4180C9", "CLUST" = "#4180C9",
    "Clustering coefficient (t.10%)"="#4180C9","CLUSTERING0.1"="#4180C9",
    "Clustering coefficient (t.20%)"="#4180C9","CLUSTERING0.2"="#4180C9", 
    "Clustering coefficient (t.30%)"="#4180C9","CLUSTERING0.3"="#4180C9",
    "Global efficiency" = "#3468A4","GLOB_EFF" = "#3468A4",
    "Global efficiency (t.10%)" = "#3468A4", "GLOB_EFF0.1" = "#3468A4",
    "Global efficiency (t.20%)" = "#3468A4", "GLOB_EFF0.2" = "#3468A4",
    "Global efficiency (t.30%)" = "#3468A4","GLOB_EFF0.3" = "#3468A4",
    "Persistent Homology" = "#95435C",  
    "Backbone Strength" = "#95435C", "PH_BS" = "#95435C",
    "Backbone Dispersion" = "#733447","PH_BD" = "#733447",
    "Cycle Strength" = "#522633","PH_CS" = "#522633", 
    "Minimum Spanning Tree" = "#C0915C",
    "Diameter" = "#DCA769","MST_DIAM" = "#DCA769",
    "Leaf fraction" = "#C0915C","MST_LEAF" = "#C0915C",
    "Raw functional connectivity" = "#9933FF",
    "Mean connectivity" = "#9933FF","RAW_FC" = "#9933FF"
  )
  
  #adapt colours depending on number first label in the x label
  df_contrasts$prefix <- sapply(df_contrasts$contrast_label, function(label) {
    matched <- grep(paste0("^", names(color_map), collapse = "|"), label, value = TRUE)
    if (length(matched) == 0) return(NA)
    matched_prefix <- names(color_map)[sapply(names(color_map), function(p) grepl(paste0("^", p), label))]
    if (length(matched_prefix) > 0) matched_prefix[1] else NA
  })
  df_contrasts$color <- color_map[df_contrasts$prefix]
  
  #categorise contrasts broadly by fMRI task
  df_contrasts$taskFC_leader <- str_match(df_contrasts$contrast, "^[^ ]+ ([^ ]+(?: \\(\\dT\\))?)")[,2]
  df_contrasts <- df_contrasts %>%
    arrange(taskFC_leader, label_leader)
  df_contrasts$contrast_label <- factor(df_contrasts$contrast_label, levels = df_contrasts$contrast_label)
  #add dashed lines between fMRI tasks
  group_breaks <- df_contrasts %>%
    group_by(taskFC_leader) %>%
    summarise(last = last(contrast_label)) %>%
    mutate(y = match(last, df_contrasts$contrast_label) + 0.5)  #offset to draw line after last contrast of that fMRI task
  #add the corresponding fMRI task label too
  taskFC_labels <- df_contrasts %>%
    mutate(row = row_number()) %>%
    group_by(taskFC_leader) %>%
    summarise(x_pos = mean(row))
  taskFC_labels$y_pos <- 0.8
  #rename fMRI task for clarity
  taskFC_truenames <- c( "WM" = "WORKING\nMEMORY", "LG" = "LANGUAGE", "SOCIAL" = "SOCIAL\nCOGNITION")
  taskFC_labels <- taskFC_labels %>%
    mutate(taskFC_leader = recode(taskFC_leader, !!!taskFC_truenames))
  
  #rename legend for clarity
  category_names <- c( "GLOB_EFF" = "Global efficiency", "CLUST" = "Clustering coefficient", "PH_BS" = "Backbone Strength", "PH_CS" = "Cycle Strength", "PH_BD" = "Backbone Dispersion", "MST_LEAF" = "Leaf fraction", "MST_DIAM" = "Diameter", "RAW_FC" = "Mean connectivity" )
  df_contrasts <- df_contrasts %>%
    mutate(prefix = recode(prefix, !!!category_names))
  #only plot at for second plot to avoid repeating
  if (scan == "7T") {
    legend_theme <- theme(
      legend.position = "bottom",
      legend.justification = c(0, 0.5),
      legend.box = "horizontal"
    )
    legend_scale <- scale_color_manual(values = color_map, name = "Stronger network measure:")
    plot_title=NULL
    y_label="Estimate ± SE"
  } else {
    legend_theme <- theme(legend.position = "none")
    legend_scale <- scale_color_manual(values = color_map, guide = "none")
    plot_title="Contrasts in estimated log(partial R-squared) between conditions"
    y_label=NULL
  }
  
  #plot
  meascontrast_plot=ggplot(df_contrasts, aes(x = fct_inorder(contrast_label), y = estimate)) +
    geom_point(aes(color = prefix), size = 3) +  
    geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE, color = prefix), width = 0.2, linewidth = 1.2) +
    geom_text(aes(label = sig, y = estimate + 1.1 * SE), # nudges asterisk above the bar
              hjust = 0, size = 6, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    coord_flip() +
    labs(y=y_label, x = "",
         title=plot_title) + 
    scale_color_manual(values = color_map, name='Stronger network measure') +
    theme_minimal() +
    theme(
      plot.title.position = "plot", 
      plot.title = element_text(hjust = 0.5)) +
    legend_theme +
    geom_vline(data = group_breaks, aes(xintercept = y), linetype = "dashed", color = "grey70", linewidth = 0.6) +
    geom_text(data = taskFC_labels,
              aes(x = x_pos, y = y_pos, label = taskFC_leader),
              inherit.aes = FALSE,
              angle = 0,
              hjust = 0.5,
              color = "grey30", alpha = 0.15,
              size = 6, fontface = "bold",
              family = "Impact") +
      guides(color = guide_legend(nrow = 1))
  
  assign(x = paste0('plot_',scan), meascontrast_plot)
}

layout <- matrix(c(1, 2), ncol = 1, byrow = TRUE)
finalplot <- gridExtra::grid.arrange(
  plot_3T, plot_7T,
  layout_matrix = layout,
  heights = c(0.862, 0.138)
)
plot(finalplot)

ggsave(filename = 'figures/taskFC/pairwise_comparisons_restbaseline.png', plot = finalplot, units = 'px', width = 4200, height = 6800, dpi = 300)