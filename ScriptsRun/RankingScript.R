library(dplyr)
library(votesys)

# Put all datasets into a list
all_datasets <- list(egeod86618_results_all[[1]],
                     egeod111727_results[[1]], thca_gse_154763_results[[1]], 
                     emtab9221_results[[1]], Echad35_results[[1]])


# Function to rank methods and split by Metrics in each dataset
rank_methods <- function(df) {
  df %>%
    group_by(Metrics, Methods) %>%
    summarise(Mean = mean(metrics)) %>%
    mutate(Rank = rank(-Mean)) %>%
    arrange(Rank) %>%
    split(.$Metrics)
}

# Rank methods in each dataset and split by Metrics
ranks <- lapply(all_datasets, rank_methods)

names(ranks)=c("egeod86618","egeod111727","thca_gse_154763","emtab9221","echad35")


####################

Metrics=c("Accuracy","Sensitivity","Specificity","F1" )
RankingResults=vector(mode = "list",length = 4)

for ( i in 1:length(Metrics)){
  #Ranking 
  raw <- c(
    ranks[[1]][[Metrics[i]]]$Methods, 
    ranks[[2]][[Metrics[i]]]$Methods, 
    ranks[[3]][[Metrics[i]]]$Methods,
    ranks[[4]][[Metrics[i]]]$Methods,
    ranks[[5]][[Metrics[i]]]$Methods
  )   
  
raw1 = unique(raw) # This gets unique values from all the rank lists
  
raw <- matrix(raw, ncol =5  , byrow = TRUE)  #  ncol is the number of datasets

vote <- create_vote(raw, xtype = 2, candidate = c(raw1))
y <- borda_method(vote,modified = TRUE,allow_dup = TRUE,min_valid = 1)

RankingResults[[i]]=as.data.frame(sort(y[["other_info"]][["count_max"]],decreasing = TRUE ) ) 
names(RankingResults)[i]=Metrics[i]

}


RankingResults
# 
#Ranking Barplot
ColorFun <- colorRampPalette( c( "#CCCCCC" , "#104E8B" ) )

for (i in seq_along(RankingResults)) {
  ColorPaleta <- ColorFun( n = nrow( x = RankingResults[[i]] ) )
  RankingResults[[i]]$Color <-
    as.character(
      x = cut(
        x = rank( x = RankingResults[[i]][,1])  # used to assign order in the event of ties
        , breaks = nrow( x = RankingResults[[i]])  # same as the 'n' supplied in ColorFun
        , labels = ColorPaleta  # label the groups with the color in ColorPaleta
      )
    )
}


for (i in seq_along(RankingResults)) {
  barplot( height = RankingResults[[i]][,1]
           , names.arg = rownames(RankingResults[[i]])
           , las = 2
           , col = RankingResults[[i]]$Color
           , border = NA  # eliminates borders around the bars
           , main = paste0("Ranking of Methods for ", names(RankingResults)[i])
           , ylab = "Count of Votes"
           , ylim = c(0,25)
           , xlab = "Methods "
           , cex.names=0.5
  )
}

####################
# 
# 
# 
# # Combine ranks from all datasets
# combined_ranks <- do.call(rbind, ranks)
# 
# # Compute Borda count for each method
# borda_count <- combined_ranks %>%
#   group_by(Metrics, Methods) %>%
#   summarise(Borda = sum(Rank)) %>%
#   mutate(Borda_Rank = rank(Borda))
# 
# 
# 
# 
# library(ggplot2)
# 
# # Sort data frame by Borda count
# borda_count_sorted <- borda_count[order(-borda_count$Borda),]
# 
library(ggplot2)
library(tidyr)
library(tibble)

# Add column names
RankingResults <- lapply(RankingResults, function(df) {
  colnames(df) <- c("Count")
  return(df)
})


# Combine all dataframes into one, adding a 'Metric' column
all_data <- do.call(rbind, lapply(seq_along(RankingResults), function(i) {
  df <- RankingResults[[i]]
  df$Metric <- names(RankingResults)[i]
  df
}))

# Convert row names to a column
all_data <- all_data %>% rownames_to_column("Methods")
colnames(all_data)[2]="Value"

# Remove numerical suffixes
all_data$Methods <- sub("[0-9]+$", "", all_data$Methods)

# Plot

ggplot(all_data, aes(x = Methods, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set2") +
  labs(x = "", y = "Count of Votes", fill = "Metrics",
       title = "Ranking of Methods") +
  theme(
    text = element_text(size=16),  # General text size, including legend
    axis.title.x = element_text(size=16, face="bold"),  # x axis title
    axis.title.y = element_text(size=16, face="bold"),  # y axis title
    axis.text.x = element_text(size=12, angle = 80, hjust = 1, vjust = 1),  # x axis text
    axis.text.y = element_text(size=12),  # y axis text
    legend.text = element_text(size=12),  # Legend text
    panel.background = element_rect(fill = "white"), 
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

