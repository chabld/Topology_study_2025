#splitter script to get 3 figures instead of 1 massive one for s7 outputs. Is run during s7 via source().
library(ggplot2)
combined_table=readRDS('datasets_stat/stat_output_HCP_alltasks.rds')
label_lookup <- combined_table %>% select(task_category, category) %>% distinct()

#Rename level for clarity
levels(combined_table$cog_name)[ levels(combined_table$cog_name) == "Self-regulation/Impulsivity (Delay Discounting 200k)" ] <- "Self-regulation (Delay Discounting 200)"
levels(combined_table$cog_name)[ levels(combined_table$cog_name) == "Processing Speed  (Pattern Completion Processing Speed)" ] <- "Processing Speed  (Pattern Comparison Processing Speed)"
levels(combined_table$cog_name)[ levels(combined_table$cog_name) == "Episodic Memory (Picture Sequence Memory)" ] <- "Visual Episodic Memory (Picture Sequence Memory)"


cog_levels <- levels(combined_table$cog_name)
cog_groups <- split(cog_levels, ceiling(seq_along(cog_levels) / 5))
for (i in seq_along(cog_groups)) {
  cog_subset <- cog_groups[[i]]
  data_subset <- combined_table %>% filter(cog_name %in% cog_subset)
  
  #recompute facet-specific y positions
  facet_max_y <- data_subset %>%
    group_by(cog_name) %>%
    summarise(y_pos = max(rsq_vals, na.rm = TRUE), .groups = "drop")
  
  taskFC_labels <- data_subset %>%
    group_by(cog_name, taskFC) %>%
    summarise(x_pos = mean(as.numeric(task_category)), .groups = "drop") %>%
    left_join(facet_max_y, by = "cog_name")
  
  #recompute taskFC breaks for this subset
  taskFC_breaks <- data_subset %>%
    group_by(taskFC) %>%
    summarise(x_break = max(as.numeric(task_category)) + 0.5, .groups = "drop") %>%
    arrange(x_break) %>%
    mutate(rank = row_number()) %>%
    filter(rank < max(rank)) %>%
    ungroup()
  
  #rebuild the plot
  p_split <- ggplot(data_subset, aes(x = task_category, y = rsq_vals, group = ntw_type)) +
    geom_line(aes(group = interaction(taskFC, ntw_type), color = ntw_type), linewidth = 0.7) +
    geom_ribbon(aes(x = as.numeric(task_category), ymin = min_rsq, ymax = rsq_vals,
                    group = interaction(taskFC, ntw_type), fill = ntw_type), alpha = 0.2) +
    geom_point(aes(fill = ntw_type, shape = ntw_type, color = ntw_type), size = 3, stroke = 0.5) +
    #Significant asterisks, FDR corrected
    geom_point(data = subset(data_subset, p_vals_fdr < 0.05),
               aes(fill = ntw_type), shape = 8, size = 3, stroke = 1, color = "black", show.legend = FALSE) +
    geom_point(data = subset(data_subset, p_vals < 0.05),
               aes(fill = ntw_type), shape = 8, size = 0.5, stroke = 1, color = "black", show.legend = FALSE) +
    geom_text(aes(label = paste0('β=', round(adjb_vals, 2)),
                  y = rsq_vals + y_offset, hjust = label_hjust),
              size = 3, color = "#5E5E5E", angle = 90, vjust = 0.5) + #beta values
    geom_vline(data = taskFC_breaks, aes(xintercept = x_break),
               linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_text(data = taskFC_labels,
              aes(x = x_pos, y = y_pos, label = taskFC),
              inherit.aes = FALSE, color = "grey30", alpha = 0.15,
              size = 6, fontface = "bold", family = "Impact") +
    scale_x_discrete(labels = setNames(label_lookup$category, label_lookup$task_category)) +
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
    labs(
      x = "Network topology measure",
      y = "Partial R-Squared",
      fill = "Category",
      shape = "Category",
      color = "Category"
    ) +
    facet_wrap(~cog_name, ncol = 1, scales = "free_y",
               labeller = label_wrap_gen(width = Inf)) +
    theme_minimal() +
    theme(
      plot.margin = margin(10, 10, 10, 20),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_text(size = 16, vjust=1),
      axis.text.x = element_text(size=12, angle = 90, hjust = 1),
      strip.text = element_text(size = 16, face = "bold", lineheight = 1.3, margin = margin(t = 6, b = 6)),
      legend.position = "bottom",
      legend.justification = "centre",
      legend.margin = margin(-10, 0, 0, 0)
    ) +
    guides(fill = guide_legend(override.aes = list(stroke = 0)))
  
  #save each plot
  ggsave(
    filename = paste0("figures/taskFC/HCP_task-specific_perf_col", i, ".png"),
    plot = p_split,
    units = "px", width = 5000, height = 6000, limitsize = FALSE
  )
}
