#This script collects all occurences of a network measure having a specified predictive performance (as per significance criterion) from the set of linear models run in s7
#Chi Square and residuals post-hoc tests are used to evaluate which network measures are more frequently the best performing predictor specifically

library(chisquare)  
library(dplyr)
library(tidyr)
library(tidyverse)

###############################################################
#######################Formatting data#########################
###############################################################

for (scan in c('3T','7T'))
{
  #load data from cog specific linear models (s7)
  overall_model=readRDS('datasets_stat/stat_output_HCP_alltasks.rds')
  
  #only applicable tasks
  if(scan=="3T")
  {
    overall_model=overall_model[which(overall_model$taskFC != 'REST (7T)' & overall_model$taskFC != 'MOVIE'),]
  } else if (scan=="7T")
  {  
    overall_model=overall_model[which(overall_model$taskFC == 'REST (7T)' | overall_model$taskFC == 'MOVIE'),]
  }
  
  #choose a criterion to classify Chisquare cells
  ## flag those that are significant at p < .05
  overall_model$is_max='N'; overall_model$is_max[which(overall_model$p_vals <= 0.05)]='Y'
  
  ####################################################################
  ##########################Prepare tables############################
  ####################################################################

  contingency_table <- overall_model %>%
    count(is_max, category) %>%
    group_by(is_max, category) %>% # Ensure proper aggregation
    summarise(n = sum(n), .groups = "drop") %>% 
    filter(!is.na(is_max)) %>% # Remove rows with missing values
    pivot_wider(names_from = category, values_from = n, values_fill = list(n = 0)) %>% # Set row names to "Y" and "N"
    column_to_rownames(var = "is_max") %>%
    as.matrix() 
  
  #exclude variables with empty cells
  for (var in colnames(contingency_table))
  { 
    if(0 %in% contingency_table[,var]) 
    {
     message(paste(var,"in", scan, "has 0 significant prediction, so cannot even be tested in the chi square, will be reported as NA"))
     contingency_table <- contingency_table[ , !(colnames(contingency_table) %in% var)]
    }
  }
  
  ####################################################################
  ####################Running Chisquare test##########################
  ####################################################################

  chisq <- chisquare(data = contingency_table, 
                     reference.level = 1, 
                     row.level = 2)
  
  #chisquare model's p-value 
  pval=chisq$chi.sq.related.results$chisq.p.value.MC
  
  print(paste0('For ',scan,' scans:'))
  print(paste0('X2=',chisq$chi.sq.related.results$chisq.statistic, '; Cohen\'s W=', chisq$chi.sq.based.assoc.measures$W, "; p=",pval ))
  
  #Get adjusted residuals
  adj_res=chisq$post.hoc$adj.stand.resid
  adj_res_p = 2 * (1 - pnorm(abs(adj_res))) #convert them to pvalues
  #store for report:
  assign(paste0('adj_res_',scan), adj_res)
  assign(paste0('adj_res_p_',scan), adj_res_p)
  
  #moment-corrected residuals
  momcor_res=chisq$post.hoc$mom.corr.stand.resid  
  momcor_res_p = 2 * (1 - pnorm(abs(momcor_res))) #convert them to pvalues
  #store for report:
  assign(paste0('momcor_res_',scan), momcor_res)
  assign(paste0('momcor_res_p_',scan), momcor_res_p)
  
}

####################################################################
####################Report and FDR corrections######################
####################################################################

#final report for both 3T and 7T scans:
all_momcor_res <- bind_rows(as.data.frame(momcor_res_3T), as.data.frame(momcor_res_7T))
#get p-values across scanners
all_momcor_res_p=all_momcor_res %>%
  mutate(across(where(is.numeric), ~ 2 * (1 - pnorm(abs(.)))))
all_momcor_res_p=all_momcor_res_p[c(2,4),]
row.names(all_momcor_res_p)=c('Y_3T','Y_7T')
#only correct across positive differences (one-tailed)
all_momcor_res_p_fdr=p.adjust(unlist(all_momcor_res_p), method = 'fdr') 
all_momcor_res_p_fdr <- as.data.frame(matrix(
  format(p.adjust(unlist(all_momcor_res_p), method = "fdr"), scientific=F),
  nrow = nrow(all_momcor_res_p),
  byrow = FALSE
))
colnames(all_momcor_res_p_fdr) <- colnames(all_momcor_res_p)
rownames(all_momcor_res_p_fdr) <- rownames(all_momcor_res_p)

print('Moment corrected standardised residuals (sig. at p(FDR) < .05):')
print(rbind(colnames(all_momcor_res_p_fdr)[col(all_momcor_res_p_fdr)[all_momcor_res_p_fdr <= 0.05]],
      all_momcor_res_p_fdr[all_momcor_res_p_fdr <= 0.05]
))

#in article, all_momcor_res was transposed with t(all_momcor_res) for readability
