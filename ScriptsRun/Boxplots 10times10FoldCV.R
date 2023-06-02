library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(tidyverse)
library(RColorBrewer)


egeod111727_results[[1]]
egeod86618_results_all[[1]]
results_all_thca_gse_154763[[1]]
Results_emtab9221[[1]]
ResultsEchad35[[1]]

#Mergelist <- bind_rows(egeod86618_results_all[[2]])

data_list <- unlist(ResultsEchad35[[2]])


# Convert the list to a data frame
accuracy_df <- data_list %>% 
  enframe(name = "Method", value = "Accuracy") %>% 
  unnest(cols = c(Accuracy))

accuracy_df <- accuracy_df %>%
  separate(Method, into = c("Method", "ML_Method"), sep = "\\.(?=[^.]*$)", remove = FALSE)


accuracy_df <- accuracy_df %>%
  mutate(ML_Method = str_remove(ML_Method, "[0-9]+"))



# Determine the number of unique methods
n_methods <- accuracy_df %>% 
  distinct(Method) %>% 
  nrow()

accuracy_df$ML_Method=toupper(accuracy_df$ML_Method)
accuracy_df$ML_Method=ifelse(accuracy_df$ML_Method=="NAIVE_BAYES","NAIVE BAYES",accuracy_df$ML_Method)

# Create a color palette with the required number of colors
my_palette <- colorRampPalette(brewer.pal(8, "Set3"))(n_methods)

BoxPlot=ggplot(accuracy_df, aes(x = Method, y = Accuracy, fill = Method)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.3, alpha = 0.3, size = 0.5) +
  scale_fill_manual(values = my_palette) +
  labs(title = "ehcad_35 ~ 10 times of 10-fold CV", 
       x = "", 
       y = "") +
  facet_wrap(~ML_Method, scales = "free", ncol = 2) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 9, face = "bold", angle = 0, hjust = 1),
    axis.text.y = element_text(size = 9, face = "bold"),
    panel.grid.major = element_line(color = "gray", linetype = "dotted"),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black", size = 1),
    strip.text = element_text(face = "bold") # Change facet label text
  ) +
  coord_flip() # To make horizontal boxplots for better visualization if you have many categories

BoxPlot
BoxPlot1=BoxPlot
BoxPlot2=BoxPlot
BoxPlot3=BoxPlot
BoxPlot4=BoxPlot
BoxPlot5=BoxPlot

# ~10 times 10-fold CV

library(cowplot)  
# Use plot_grid to arrange the plots
final_plot <- plot_grid(BoxPlot1, BoxPlot2, BoxPlot3, BoxPlot4, BoxPlot5, 
                        ncol = 2, nrow = 3, 
                       # align = "h",# this aligns the plots vertically
                        # this adjusts the height of each plot
                        rel_widths = 1,
                        rel_heights = c(4, 4, 4),
                        labels = "AUTO") 
print(final_plot)

##### OTher Boxplot all to one 



egeod111727_results[[1]]
egeod86618_results_all[[1]]
results_all_thca_gse_154763[[1]]
Results_emtab9221[[1]]
ResultsEchad35[[1]]

#Mergelist <- bind_rows(egeod86618_results_all[[2]])

data_list <- unlist(egeod111727_results[[2]])


# Convert the list to a data frame
accuracy_df <- data_list %>% 
  enframe(name = "Method", value = "Accuracy") %>% 
  unnest(cols = c(Accuracy))

accuracy_df <- accuracy_df %>%
  separate(Method, into = c("Method", "ML_Method"), sep = "\\.(?=[^.]*$)", remove = FALSE)


accuracy_df <- accuracy_df %>%
  mutate(ML_Method = str_remove(ML_Method, "[0-9]+"))



# Determine the number of unique methods
n_methods <- accuracy_df %>% 
  distinct(Method) %>% 
  nrow()

accuracy_df$ML_Method=toupper(accuracy_df$ML_Method)
accuracy_df$ML_Method=ifelse(accuracy_df$ML_Method=="NAIVE_BAYES","NAIVE BAYES",accuracy_df$ML_Method)


# Create a color palette with the required number of colors
my_palette <- colorRampPalette(brewer.pal(12, "Set3"))(n_methods)

BoxPlot <- ggplot(accuracy_df, aes(x = Method, y = Accuracy, fill = ML_Method)) +
  geom_boxplot() +
  #geom_jitter(alpha = 0.2,size = 1,  position = position_jitter(width = 0.2)) +
  labs(x = " ", y = "10 times of 10 folds CV") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme( 
         axis.text.x =element_text(angle = 60, hjust = 1), # axis.text.x = element_blank(), # 
        strip.background = element_rect(fill = "white", color = "black", size = 1),
        strip.text = element_text(face = "bold"),
        panel.border = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"
       ) )


BoxPlot1=BoxPlot
BoxPlot2=BoxPlot
BoxPlot3=BoxPlot
BoxPlot4=BoxPlot
BoxPlot5=BoxPlot



library(cowplot)  
# Use plot_grid to arrange the plots
final_plot <- plot_grid(BoxPlot1, BoxPlot2, BoxPlot3, BoxPlot4, BoxPlot5, 
                        ncol = 2, nrow = 3, 
                         align = "hv",# this aligns the plots vertically
                        # this adjusts the height of each plot
                        #rel_widths = 1,
                        #rel_heights = c(1, 2, 2),
                        scale = 0.9,
                        labels = "AUTO") 
print(final_plot)


