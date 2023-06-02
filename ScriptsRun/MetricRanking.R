# Assuming that the 'metrics' column in your data is the performance measure (higher is better)
# And that 'Methods' contains the feature selection methods

egeod86618_results_all[[1]]$rank <- with(egeod86618_results_all[[1]], ave(metrics, Methods, FUN = rank))
egeod111727_results[[1]]$rank <- with(egeod111727_results[[1]], ave(metrics, Methods, FUN = rank))
thca_gse_154763_results[[1]]$rank <- with(thca_gse_154763_results[[1]], ave(metrics, Methods, FUN = rank))
emtab9221_results[[1]]$rank <- with(emtab9221_results[[1]], ave(metrics, Methods, FUN = rank))
Echad35_results[[1]]$rank <- with(Echad35_results[[1]], ave(metrics, Methods, FUN = rank))




ggplot(ResultsEchad34[[1]], aes(x = Methods, y = -rank, color = ML_Method)) +
  geom_point() +
  facet_wrap(~Metrics) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Ranking of Feature Selection Methods", 
       x = "Feature Selection Method", 
       y = "Rank")



# First add a new column in each data frame indicating the dataset
egeod86618_results_all[[1]]$dataset <- "egeod86618"
egeod111727_results[[1]]$dataset <- "egeod8661"
thca_gse_154763_results[[1]]$dataset <- "thca_gse_154763"
emtab9221_results[[1]]$dataset <- "emtab9221"
Echad35_results[[1]]$dataset <- "echad35"


# Combine the datasets
DFranking <- rbind(egeod86618_results_all[[1]], egeod111727_results[[1]], 
                   thca_gse_154763_results[[1]], emtab9221_results[[1]], Echad35_results[[1]])


a=DFranking %>%
  group_by(Metrics, ML_Method,dataset) %>%
  summarise(Mean = mean(metrics)) %>%
  split(.$Metrics)

a


ggplot(DFranking, aes(x = Methods, y = rank, color = ML_Method)) +
  geom_point() +
  facet_wrap(~Metrics + dataset) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Ranking of Feature Selection Methods", 
       x = "Feature Selection Method", 
       y = "Rank")




ggplot(DFranking, aes(x = Methods, y = rank, color = ML_Method)) +
  geom_point(size = 1) +
  facet_wrap(~Metrics + dataset, strip.position = "bottom") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Ranking of Feature Selection Methods across Datasets", 
    x = "Feature Selection Method", 
    y = "Rank"
  )
DFranking



# Ranking the methods by metric
d=DFranking %>%
  group_by(Metrics) %>%
  mutate(rank = rank(metrics     )) %>%
  arrange(Metrics, rank)



ggplot(d, aes(x = Methods, y = rank, color = ML_Method)) +
  geom_point(size = 1) +
  facet_wrap(~Metrics + dataset, strip.position = "bottom") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Ranking of Feature Selection Methods across Datasets", 
    x = "Feature Selection Method", 
    y = "Rank"
  )


# Calculate the mean of each metric for each method
mean_metrics <- aggregate(ResultsEchad34[[1]]$metrics ~ ResultsEchad34[[1]]$Methods + ResultsEchad34[[1]]$Metrics, FUN = mean)

# Rename columns for better readability
colnames(mean_metrics) <- c("Methods", "Metrics", "Mean")

# Order the data frame by the mean of each metric
mean_metrics <- mean_metrics[order(mean_metrics$Metrics, -mean_metrics$Mean),]

# Display the top method for each metric
top_methods <- aggregate(Mean ~ Metrics, mean_metrics, function(x) { mean_metrics$Methods[which.max(x)] })



# Calculate the mean of each metric for each method
mean_metrics <- aggregate(ResultsEchad34[[1]]$metrics ~ ResultsEchad34[[1]]$Methods + ResultsEchad34[[1]]$Metrics, FUN = mean)

# Rename columns for better readability
colnames(mean_metrics) <- c("Methods", "Metrics", "Mean")

# Split the data frame by metric
metrics_list <- split(mean_metrics, mean_metrics$Metrics)

# Sort each data frame by the mean metric values
metrics_list <- lapply(metrics_list, function(df) df[order(-df$Mean), ])

# Display the list of data frames
metrics_list



###########
i <- seq(1:20)
updated_list <- lapply(egeod111727_results[[2]], function(x) {
  x[[4]] <- ifelse(x[[4]] == 1, 0.99, x[[4]])
  x
})

egeod111727_results[[2]]=updated_list
