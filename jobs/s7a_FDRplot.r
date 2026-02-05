#script taylored specifically for the results found in S7, to report only panels with at least one prediction that survives FDR correction. Will likely not adapt dynamically well to other results.

###############################################################
###############################################################
##############################PLOTTING#########################

library(ggplot2)
combined_table=readRDS('datasets_stat/stat_output_HCP_alltasks.rds') 

#FDR correction
combined_table$p_vals_fdr=p.adjust(combined_table$p_vals,method = 'fdr')
combined_table <- combined_table %>%
  group_by(taskFC, cog_name) %>%
  filter(any(p_vals_fdr < 0.05)) %>%
  ungroup()

#rename cognitive measures for clarity in the plot
cog_names <- c(
  "Fluid Intelligence (Penn Progressive Matrices)"= "Penn Progressive Matrices",    
  "Vocabulary Comprehension  (Picture Vocabulary)"="Vocabulary",         
  "Spatial Orientation (Variable Short Penn Line Orientation Test)"="Spatial\n\n\n\nOrientation",
  "Cognition Fluid Composite"="Cognition Fluid\n\n\n\n(Composite)", 
  "Cognition Crystallized Composite"="Cognition\n\n\n\nCrystallized\n\n\n\n(Composite)",
  "Verbal Episodic Memory (Penn Word Memory Test)"="Episodic Memory\n\n\n\n(Verbal)",
  "Episodic Memory (Picture Sequence Memory)"="Episodic Memory\n\n\n\n(Visual)",
  "Self-regulation/Impulsivity (Delay Discounting 200k)"="Self-regulation"
)
#apply
combined_table$cog_name <- dplyr::coalesce(
  cog_names[as.character(combined_table$cog_name)],
  as.character(combined_table$cog_name)
)

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
#add categories per taskFC 
combined_table$task_category <- paste0(combined_table$category, "_", combined_table$taskFC)
#create lookup table for category ordering within ntw_type
category_order <- tibble::tibble(
  category_abbr = names(category_names),
  category = unname(category_names),
  ntw_type = dplyr::case_when(
    category_abbr %in% c("CLUSTERING0.1", "CLUSTERING0.2", "CLUSTERING0.3", "GLOB_EFF0.1", "GLOB_EFF0.2", "GLOB_EFF0.3") ~ "Graph measures",
    category_abbr %in% c("MST_DIAM", "MST_LEAF") ~ "Minimum Spanning Tree",
    category_abbr %in% c("PH_BD", "PH_BS", "PH_CS") ~ "Persistent Homology",
    category_abbr == "RAW_FC" ~ "Raw functional connectivity"
  )
) %>%
  group_by(ntw_type) %>%
  arrange(category) %>%  #alphabetically ordered in each ntw_type
  mutate(cat_order = row_number()) %>%
  ungroup()

#reorder task_category accordingly
combined_table <- combined_table %>%
  left_join(category_order, by = c("category", "ntw_type")) %>%
  arrange(taskFC, ntw_type, cat_order) %>%
  mutate(task_category = factor(task_category, levels = unique(task_category)))

#get x-axis positions for each fMRI task block
taskFC_breaks <- combined_table %>%
  group_by(cog_name, taskFC) %>%
  summarise(x_break = max(as.numeric(task_category)) + 0.5, .groups = "drop") %>%
  group_by(cog_name) %>%
  mutate(rank = row_number()) %>%
  filter(rank < max(rank)) %>%
  ungroup()

#this adds the offset for relative position of beta value
combined_table <- combined_table %>%
  group_by(cog_name) %>%
  mutate(y_offset = mean(rsq_vals) + 0.001)  
combined_table <- combined_table %>% mutate(label_hjust = 0.5) 

#dropped non showing levels after FDR correction
combined_table <- combined_table %>%
  droplevels()

#ensures each facet has its own local x-axis scale:
combined_table <- combined_table %>%
  group_by(cog_name) %>%
  mutate(task_category = factor(task_category, levels = unique(task_category))) %>%
  ungroup()

#main plot
p3 <- ggplot(combined_table, aes(x = task_category, y = rsq_vals, group = ntw_type)) +
  geom_line(aes(group = interaction(taskFC, ntw_type), color = ntw_type), linewidth = 0.7) +
geom_ribbon(aes(x = task_category, 
                ymin = min_rsq, ymax = rsq_vals, 
                group = interaction(taskFC, ntw_type), fill = ntw_type), 
            alpha = 0.2) + 
  geom_point(aes(fill = ntw_type, shape = ntw_type, color = ntw_type), size = 3, stroke = 0.5) +
  #Significant asterisks, FDR corrected
  geom_point(data = subset(combined_table, p_vals_fdr < 0.05), 
             aes(fill = ntw_type), shape = 8, size = 3, stroke = 1, color = "black", show.legend = FALSE) + 
  #not FDR corrected: lower size
  geom_point(data = subset(combined_table, p_vals < 0.05), 
             aes(fill = ntw_type), shape = 8, size = 0.5, stroke = 1, 
             color = "black", show.legend = FALSE) +
  #adjusted beta values
  geom_text(aes(label = paste0('β=', round(adjb_vals, 2)),
                y = rsq_vals + y_offset,
                hjust = label_hjust),
            size = 4, color = "#5E5E5E",
            angle=90, vjust=0.5) +
  labs(
    x = "Network Measure",
    y = "Partial R-Squared",
    #title = paste0("Predicting ", response_var, " from network topology measures in resting-state fMRI"),
    fill = "Category",
    shape = "Category",
    color = "Category"
  ) +
  facet_wrap(~cog_name, nrow = 1, scales = "fixed", labeller = label_wrap_gen(width = Inf))+
  facet_grid(~ cog_name, scales = "free_x", space = "free_x")+
  theme_minimal() +
  theme(strip.text = element_text(
    size = 10,
    face = "bold",
    lineheight = 0.2,
    margin = margin(t = 6, b = 6)  #top and bottom padding
  )) +
  theme(plot.margin = margin(10, 10, 10, 20),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 16, vjust=1),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=12, angle = 90, hjust = 1),
        strip.text = element_text(size = 12, face = "bold", 
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
    legend.margin = margin(-10, 0, 0, 0),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14)) +
  guides(fill = guide_legend(override.aes = list(stroke = 0)))  #legend clean, no halos

multi_taskFC_facets <- combined_table %>%
  group_by(cog_name) %>%
  summarise(n_taskFC = n_distinct(taskFC), .groups = "drop") %>%
  filter(n_taskFC > 1) %>%
  pull(cog_name)
multi_panel <- combined_table %>% filter(cog_name %in% multi_taskFC_facets)

#add dashed lines between taskFC
p3 <- p3 +
  geom_vline(data = taskFC_breaks,
             aes(xintercept = x_break),
             linetype = "dashed", color = "grey50", linewidth = 0.5)

#print taskFC name between the dashed lines
y_fixed <- max(combined_table$rsq_vals, na.rm = TRUE) - 0.01  #offset
global_y_mid <- max(combined_table$rsq_vals, na.rm = TRUE) / 2
taskFC_labels <- combined_table %>%
  group_by(cog_name, taskFC) %>%
  summarise(
    task_category = task_category[which.min(abs(as.numeric(task_category) - mean(as.numeric(task_category))) )],
    .groups = "drop"
  ) %>%
  mutate(y_pos = global_y_mid)
p3 <- p3 +
  geom_text(data = taskFC_labels,
            aes(x = task_category, y = y_pos, label = taskFC),
            inherit.aes = FALSE,
            color = "grey30", alpha = 0.15,
            size = 6, fontface = "bold",
            family = "Impact")

#don't need tasKFC name in the plot's x labels as taskFC named in the facet
label_lookup <- combined_table %>%
  group_by(task_category) %>%
  summarise(category = first(category), .groups = "drop")
p3 <- p3 + scale_x_discrete(labels = setNames(label_lookup$category, label_lookup$task_category))
print(p3)

ggsave(paste0("figures/taskFC/HCP_task-specific_perf_alltasks_FDRonly.png"), 
       plot = p3, units = 'px', width=6500, height=2400, limitsize=F)
