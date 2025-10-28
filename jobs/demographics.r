#demographics

library(dplyr)

######################################################################
######standardising subjects across fMRI task and neuropsych task#####
######################################################################

#identify subjects that have data for all task FC in 3T
tasklist=c('REST3T','MOTOR','GAMBLING','LG','WM','SOCIAL','RELATIONAL')
for (taskFC in tasklist)
{
  dataset_lm_resp=readRDS(paste0('datasets_model_data/HCPmodel_data_',taskFC,'.rds'))
  #assign
  assign(paste0('subs_',taskFC), as.character(unique(dataset_lm_resp$participant_id)))
}
common_subs <- Reduce(intersect, list(subs_GAMBLING, subs_LG, subs_MOTOR, subs_RELATIONAL, subs_REST3T, subs_SOCIAL, subs_WM))
#identify subject that have data for all neuropsych measures
toexclude=c()
for (i in 1:nrow(dataset_lm_resp))
{ if(all(!is.na(dataset_lm_resp[i,]))==FALSE)
{toexclude=c(toexclude,as.character(dataset_lm_resp$participant_id[i]))}
}
common_subs3T=common_subs[which(! common_subs %in% toexclude)]

#same of 7T scanners
tasklist=c('REST7T','MOVIE')
for (taskFC in tasklist)
{
  dataset_lm_resp=readRDS(paste0('datasets_model_data/HCPmodel_data_',taskFC,'.rds'))
  #assign
  assign(paste0('subs_',taskFC), as.character(unique(dataset_lm_resp$participant_id)))
}
common_subs <- Reduce(intersect, list(subs_REST7T, subs_MOVIE))
#identify subject that have data for all neuropsych measures
toexclude=c()
for (i in 1:nrow(dataset_lm_resp))
{ if(all(!is.na(dataset_lm_resp[i,]))==FALSE)
{toexclude=c(toexclude,as.character(dataset_lm_resp$participant_id[i]))}
}
common_subs7T=common_subs[which(! common_subs %in% toexclude)]

#force all subs to be cross-tasks
beh_data <- read.csv('HCP/HCP_behdata.csv')
beh_data3T=beh_data[which(beh_data$Subject %in% common_subs3T),]
beh_data7T=beh_data[which(beh_data$Subject %in% common_subs7T),]

library(ggplot2)
library(dplyr)
library(forcats)
library(scales)

# Combine both datasets with a scanner label
beh_data3T$scanner <- "3T scanner"
beh_data7T$scanner <- "7T scanner"

combined_data <- bind_rows(beh_data3T, beh_data7T)
# Ensure age bins are ordered consistently
age_levels = c('22-25','26-30','31-35','36+')
combined_data$Age <- factor(combined_data$Age, levels = age_levels, ordered = TRUE)

# Plot
demogplot=ggplot(combined_data, aes(x = Age, fill = Gender)) +
  geom_bar(position = "dodge") +
  facet_wrap(~scanner, ncol = 2) +
  scale_fill_manual(values = c("F" = "#009E73", "M" = "#56B4E9")) +  # colorblind-friendly
  labs(x = "Age bins", y = "Count", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(filename = 'figures/demographics.png', plot = demogplot,
       units = 'px', height=1500, width=3000, dpi = 300)

##comparisons
age_table <- table(combined_data$Age, combined_data$scanner)
chisq.test(age_table)
sex_table <- table(combined_data$Gender, combined_data$scanner)
chisq.test(sex_table)


