
library(gridExtra)
library(ggplot2)
library(tidyr)
####Heatmap Function


#### Two to One OR separate #### 


ComparisonPlotsFun <- function(Results_Conf_all) {
  
  # Convert the data frame into wide format
  results_wide <- tidyr::spread(Results_Conf_all, key = ML_Method, value = metrics)
  
  # Rearrange the columns for better visualization
  results_wide <- results_wide[, c("Methods", "Metrics", "knn", "naive_bayes", "svmLinear", "rf")]
  
  # Convert the data frame into long format
  results_heatmap <<- reshape2::melt(results_wide, id.vars = c("Methods", "Metrics"), variable.name = "ML_Method", value.name = "Value")
  results_heatmap$ML_Method=toupper(results_heatmap$ML_Method)
  results_heatmap$ML_Method=ifelse(results_heatmap$ML_Method=="NAIVE_BAYES","NAIVE BAYES",results_heatmap$ML_Method)
  plot1 <- ggplot(results_heatmap, aes(x = ML_Method, y = Metrics, fill = Value)) +
    geom_tile(color = "black", size = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Value)), color = "black", size = 3) +
    facet_wrap(~ Methods) +
    scale_fill_gradient(low = "white", high = "red") +
    labs(
         x = "",
        y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "white", color = "black", size = 1),
          strip.text = element_text(face = "bold"))
  
  plot2 <- ggplot(results_heatmap, aes(x = Methods, y = Value, color = Metrics, group = Metrics)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ ML_Method) +
    #labs(#title = "Performance Metrics per Methods and ML Method",
    #     x = "Methods",
     #    y = "Value") +
    theme_minimal() +
    scale_color_brewer(palette = "Set1") +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) # Set x-axis labels to be vertical
  
  return(list(plot1 = plot1, plot2 = plot2))
}

#egeod86618_results_all  MegeListENd   results_all   results_all_thca_gse_154763   Results_emtab9221
# egeod86618_results_all[[1]]$ML_Method=tolower(egeod86618_results_all[[1]]$ML_Method)
# results_alleEgeod_8661[[1]]$ML_Method=tolower(results_alleEgeod_8661[[1]]$ML_Method)
# results_all_thca_gse_154763[[1]]$ML_Method=tolower(results_all_thca_gse_154763[[1]]$ML_Method)
# Results_emtab9221[[1]]$ML_Method=tolower(Results_emtab9221[[1]]$ML_Method)
# MegeListENd[[1]]$ML_Method=tolower(MegeListENd[[1]]$ML_Method)


plots1 <- ComparisonPlotsFun(egeod86618_results_all[[1]])
plots2 <- ComparisonPlotsFun(results_alleEgeod_8661[[1]])
plots3 <- ComparisonPlotsFun(results_all_thca_gse_154763[[1]])
plots4 <- ComparisonPlotsFun(Results_emtab9221[[1]])
plots5 <- ComparisonPlotsFun(ResultsEchad34[[1]])

plots1$plot1
plots2$plot1
plots3$plot1
plots4$plot1
plots5$plot1


grid.arrange(plots1$plot1, plots2$plot1,plots3$plot1,plots4$plot1,plots5$plot1, ncol = 2,nrow=3,top = "Your Title")

library(grid)
library(gridExtra)

# Create an empty plot
empty <- ggplot() + theme_void()

# Arrange the plots in a custom layout
grid.newpage()
grid.draw(cbind(
  arrangeGrob(plots1$plot1, plots2$plot1, plots3$plot1, nrow = 3),
  arrangeGrob(plots4$plot1, plots5$plot1, empty, nrow = 3),
  size = "last"
))




names(results_all_thca_gse_154763[[2]])=methods

results_all[[2]]=results_all[[2]][1:3]
# Convert the list into a data frame
# Load the ggplot2 package
library(dplyr)

results_all[[2]]

# Convert the list into a data frame
accuracy_list <- Results_emtab9221[1]


# Convert the accuracy list to a long-format dataframe
accuracy_df <- bind_rows(accuracy_list, .id = "Method")

# Gather the ML_Method and Accuracy columns into a long-format dataframe
accuracy_long <- accuracy_df %>%
  tidyr::gather(key = "ML_Method", value = "Accuracy", -Method)
accuracy_long <-na.omit(accuracy_long)

library(ggplot2)

ggplot(accuracy_long, aes(x = Method, y = Accuracy, fill = ML_Method)) +
  theme_minimal() +
  # geom_jitter(position = position_jitterdodge(dodge.width = 0.8), alpha = 0.5, size = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) + # Set a smaller width value
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11)
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "10 times 10-fold CVs",
       x = "ML Method",
       y = "Accuracy") +
  scale_x_discrete(expand = c(0.1, 0))



# Aggregate the results and calculate the rankings
aggregated_data <- accuracy_long %>%
  group_by(Method, ML_Method) %>%
  summarize(MeanAccuracy = mean(Accuracy)) %>%
  mutate(Rank = rank(-MeanAccuracy)) %>%
  arrange(Method, Rank)

# Create the plot
ggplot(aggregated_data, aes(x = interaction(Method, ML_Method, sep = "_"), y = MeanAccuracy, fill = ML_Method)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = Rank), vjust = -0.5) +
  labs(x = "Method & ML Method", y = "Mean Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_blank()) +
  scale_fill_brewer(palette = "Set2")





# Create the box plot
# Create the box plot with jitter points
ggplot(accuracy_df, aes(x = Method, y = Accuracy)) +
  geom_boxplot(fill = "dodgerblue", outlier.color = "red", alpha = 0.7, color = "black", lwd = 0.6) +
  geom_jitter(aes(color = Method), width = 0.15, size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Accuracy by Method", x = "Method", y = "Accuracy") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.position = "none") +
  scale_color_manual(values = c("darkorange", "forestgreen", "purple"))

library(ggpubr)
ggboxplot(accuracy_df,x="Method",y="Accuracy", width=0.6, add="jitter",ylim = c(0.5,1), 
          title = " 10 times 10-fold CVs " ,
          color = "Method",palette = c("#5D2B00","#004A25","#8F0000","#55008F")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.5),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))


