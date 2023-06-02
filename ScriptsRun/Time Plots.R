# Please install the necessary packages before running the code.
# You can do this by uncommenting the following lines and running them:
# install.packages("ggplot2")
# install.packages("dplyr")

library(ggplot2)
library(dplyr)

# Assuming your data is in a dataframe called df
df <- AllFS_Time1

# 1. Compare the performance of the different methods across the datasets
df %>%
  ggplot(aes(x = Method, y = Time, fill = Method)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~ Dataset, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") +
  labs(title = "Time Comparison of Different Methods Across Datasets", 
       x = "Time (Log Scale)", 
       y = "Time")




df %>%
  ggplot(aes(x = Method, y = Time, fill = Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Time Comparison of Different Methods Across Datasets", 
       x = "Methods", 
       y = "Time",
       fill = "Method")


# 2. Explore the relationship between the number of cells in a dataset and the time taken for dimensionality reduction
df %>%
  ggplot(aes(x = Cells, y = Time)) +
  geom_point(aes(color = Method)) +
  labs(title = "Time vs Number of Cells", 
       x = "Number of Cells", 
       y = "Time",
       color = "Method")

# 3. Look at how the different methods perform relative to each other for each dataset
# First, we need to calculate the ratio of fastest to slowest time for each dataset
df_ratio <- df %>%
  group_by(Dataset) %>%
  summarize(ratio = min(Time) / max(Time))
df_ratio[1,1]="echad35 (22327 Cells)"
df_ratio[2,1]= "egeod111727 (320 Cells)"  
df_ratio[3,1]= "egeod86618 (537 Cells)"  
df_ratio[4,1]= "emtab9221 (6807 Cells)"  
# Now plot the ratios
df_ratio %>%
  ggplot(aes(x = Dataset, y = ratio, fill = Dataset)) +
  geom_col(width = 0.8) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none",  # Remove the legend
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray", linetype = "dotted"),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.caption = element_text(size = 10)
  ) +
  labs(
    title = "Ratio of Fastest to Slowest Time for Each Dataset",
    x = "Dataset",
    y = "Ratio"
  )



library(ggplot2)

ggplot(df, aes(x = Cells, y = Time, group = Method, color = Method)) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  labs(x = "Cells", 
       y = "Time (Log Scale)", 
       color = "Method",
       title = "Time of Different Methods by Cell Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "gray", linetype = "dotted"),
        legend.position = "right")

ggplot(df, aes(x = Cells, y = Time)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Method, scales = "free_y") +
  labs(x = "Sample Size (Cells)",
       y = "Execution Time (s)",
       title = "Time of Different Methods by Cell Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis labels vertically
    axis.text.y = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "gray", linetype = "dotted"),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black", size = 1),
    strip.text = element_text(face = "bold")
  )



# Convert the Method column to a factor to ensure it's treated as a categorical variable
AllFS_Time1$Method <- as.factor(AllFS_Time1$Method)

# Create a bar plot
ggplot(AllFS_Time1, aes(x = Method, y = log(Time), fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Method", y = "Execution Time", fill = "Dataset") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels to prevent overlap


ggplot(AllFS_Time1, aes(x =log(Time) , y = Method, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Execution Time (log)", y = "", fill = "Dataset") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background = element_rect(fill = "white"), # Set the background color here
        plot.title = element_text(hjust = 0.5, face = "bold")) 

install.packages("ggh4x")
library(ggh4x)


ggplot(AllFS_Time1, aes(x = Method, y = log(Time), group = Dataset, color = Dataset)) +
  geom_line(aes(linetype = Dataset), size = 1) +
  geom_point(size = 4) +  # Change size to a bigger value for larger points
  scale_color_brewer(palette = "Set2") +
  labs(x = "Feature Selection Methods", 
       y = "Execution Time (log scale)", 
       color = "Dataset") +
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
