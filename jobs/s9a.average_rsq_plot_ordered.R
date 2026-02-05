library(dplyr)
library(ggplot2)
library(stringr)

# --- keep your color_map ---
color_map <- c(
  "Graph_measures" = "#3468A4",
  "Clustering_coefficient" = "#4180C9", "CLUST" = "#4180C9",
  "Clustering coefficient (t.10%)"="#4180C9","CLUSTERING0.1"="#4180C9",
  "Clustering coefficient (t.20%)"="#4180C9","CLUSTERING0.2"="#4180C9", 
  "Clustering coefficient (t.30%)"="#4180C9","CLUSTERING0.3"="#4180C9",
  "Global_efficiency" = "#3468A4","GLOB_EFF" = "#3468A4",
  "Global efficiency (t.10%)" = "#3468A4", "GLOB_EFF0.1" = "#3468A4",
  "Global efficiency (t.20%)" = "#3468A4", "GLOB_EFF0.2" = "#3468A4",
  "Global efficiency (t.30%)" = "#3468A4","GLOB_EFF0.3" = "#3468A4",
  "Persistent_Homology" = "#95435C",  
  "Backbone Strength" = "#95435C", "PH_BS" = "#95435C",
  "Backbone Dispersion" = "#733447","PH_BD" = "#733447",
  "Cycle Strength" = "#522633","PH_CS" = "#522633", 
  "Minimum_Spanning_Tree" = "#C0915C",
  "Diameter" = "#DCA769","MST_DIAM" = "#DCA769",
  "Leaf fraction" = "#C0915C","MST_LEAF" = "#C0915C",
  "Raw_functional_connectivity" = "#9933FF",
  "Mean connectivity" = "#9933FF","RAW_FC" = "#9933FF"
)

for (scan in c("3T", "7T"))
{
  dataset=get(x = paste0("dataset_stat_",scan))
  dataset <- dataset %>% mutate(short_label = paste0(category_abbr, "_", taskFC))
  
  # --- summarise with mean + SE ---
  avg_vals <- dataset %>%
    group_by(short_label) %>%
    summarise(
      mean_val = mean(abs(rsq_vals), na.rm = TRUE),
      se_val   = sd(abs(rsq_vals), na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(color = {
      matched <- names(color_map)[str_detect(short_label, names(color_map))]
      if(length(matched) > 0) color_map[matched[1]] else "#999999"
    }) %>%
    ungroup() %>%
    arrange(desc(mean_val)) %>%
    mutate(short_label = factor(short_label, levels = short_label))
  
  # --- bar plot ordered by decreasing mean ---
  ggplot(avg_vals, aes(x = short_label, y = mean_val, fill = color)) +
    geom_col(show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                  width = 0.2, colour = "black") +
    scale_fill_identity() +
    labs(y = "Mean partial R-squared", x = NULL) +
    theme_minimal() +
    theme(
      axis.title.y = element_text(size = 20, vjust = 1),
      axis.text.y  = element_text(size = 12),
      axis.text.x  = element_text(size = 14, angle = 90, hjust = 1)
    )

  # Build legend_df only from categories present in dataset
  #rename legend for clarity
  category_names <- c( "GLOB_EFF0.1" = "Global efficiency", "GLOB_EFF0.3" = "Global efficiency", "GLOB_EFF0.2" = "Global efficiency", "CLUSTERING0.1" = "Clustering coefficient","CLUSTERING0.2" = "Clustering coefficient","CLUSTERING0.3" = "Clustering coefficient", "PH_BS" = "Backbone Strength", "PH_CS" = "Cycle Strength", "PH_BD" = "Backbone Dispersion", "MST_LEAF" = "Leaf fraction", "MST_DIAM" = "Diameter", "RAW_FC" = "Mean connectivity")
  dataset <- dataset %>% mutate(category = dplyr::recode(category, !!!category_names))
  #build
  legend_df <- data.frame(
    measure = factor(unique(dataset$category),
                     levels = names(color_map)))
  legend_df$measure=droplevels(legend_df$measure)
  
  #no legend needed for top plot
  if (scan == "7T") {
    legend_theme <- theme(legend.title= element_text(size = 14),
                          legend.text = element_text(size = 14),
                          legend.position = "bottom",
                          legend.box = "horizontal")
  } else {
    legend_theme <- theme(legend.position = "none")
  }
  
  # Main violin plot
  violinplot <- ggplot(dataset, aes(x = short_label, y = abs(rsq_vals))) +
    geom_violin(aes(fill = category),
                scale = "width", trim = FALSE,
                color = NA, width = 1.2, show.legend = FALSE) +
    geom_segment(data = avg_vals,
                 aes(x = as.numeric(short_label) - 0.25,
                     xend = as.numeric(short_label) + 0.25,
                     y = mean_val, yend = mean_val),
                 inherit.aes = FALSE,
                 colour = "black", linewidth = 1.2) +
    geom_errorbar(data = avg_vals,
                  aes(x = short_label,
                      ymin = mean_val - se_val,
                      ymax = mean_val + se_val),
                  inherit.aes = FALSE,
                  width = 0.2, colour = "black") +
    
    # Invisible points to trigger legend
    geom_point(data = legend_df,
               aes(x = 1, y = 0, fill = measure),
               shape = 21, alpha = 0) +
    
    scale_fill_manual(values = color_map, name = "Category") +
    labs(y = paste("pR2 across", scan,"models"), x = NULL) +
    theme_minimal() +
    legend_theme + 
    theme(
      plot.margin = margin(10,0,0,0),
      legend.margin = margin(-10, 0, 0, 0),
      axis.title.y = element_text(size = 20),
      axis.text.y  = element_text(size = 12),
      axis.text.x  = element_text(size = 14, angle = 45, hjust = 1)
    ) +
    guides(fill = guide_legend(
      nrow = 1,
      override.aes = list(alpha = 1, size=5, shape = 21, colour = NA) 
    )) +
    scale_y_continuous(limits = c(0, NA), expand = c(0,0))
  
  ##print(violinplot)
  #ggsave("figures/taskFC/violinplot_ordered.png",
  #plot = violinplot, units = 'px', width = 6500, height = 2400, limitsize = FALSE)
  
  assign(x = paste0('plot_',scan), violinplot)
}


library(patchwork) 
finalplot <- plot_3T / plot_7T 
plot(finalplot)
ggsave("figures/taskFC/violinplot_ordered.png", plot= finalplot, units = 'px', width = 6500, height = 3000, limitsize = FALSE)
