library(ggplot2)
library(gridExtra)


egeod86618_results_all[[1]]
egeod111727_results[[1]]
thca_gse_154763_results[[1]]
emtab9221_results[[1]]
Echad35_results[[1]]




Results_Conf_all=Echad35_results[[1]]


# Convert the data frame into wide format
results_wide <- tidyr::spread(Results_Conf_all, key = ML_Method, value = metrics)

# Rearrange the columns for better visualization
results_wide <- results_wide[, c("Methods", "Metrics", "knn", "naive_bayes", "svmLinear", "rf")]

# Convert the data frame into long format
results_heatmap <- reshape2::melt(results_wide, id.vars = c("Methods", "Metrics"), variable.name = "ML_Method", value.name = "Value")
results_heatmap$ML_Method=toupper(results_heatmap$ML_Method)
results_heatmap$ML_Method=ifelse(results_heatmap$ML_Method=="NAIVE_BAYES","NAIVE BAYES",results_heatmap$ML_Method)
plot <- ggplot(results_heatmap, aes(x = ML_Method, y = Metrics, fill = Value)) +
  geom_tile(color = "black", size = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Value)), color = "black", size = 4) +  # Increase the size for better readability
  facet_wrap(~ Methods) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    x = "",
    y = "") +
  theme_minimal() +
  theme(
    axis.text.x =   element_text(size = 14, angle = 80, hjust = 1), #element_blank(), #
    axis.text.y = element_text(size = 14),  # Increase the size for better readability
    strip.background = element_rect(fill = "white", color = "black", size = 1),
    strip.text = element_text(face = "bold", size = 16),  # Increase the size for better readability
    legend.title = element_text(size = 14),  # Increase the size for better readability
    legend.text = element_text(size = 12),  # Increase the size for better readability
    plot.title = element_text(size = 18),  # Increase the size for better readability
    plot.margin = margin(1, 1, 1, 1, "cm")  # Added some margin for better layout
  )




plot1=plot 
plot2=plot
plot3=plot
plot4=plot 
plot5=plot 


library(cowplot)  
# Use plot_grid to arrange the plots
final_plot <- plot_grid(plot1, plot2, plot3, plot4, plot5, 
                        ncol = 2, nrow = 3, 
                        #align = "h",# this aligns the plots vertically
                       # this adjusts the height of each plot
                        rel_widths = 1,
                        rel_heights = c(1, 1, 1.2),
                       labels = "AUTO") 
print(final_plot)
?plot_grid
