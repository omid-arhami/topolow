# Load necessary libraries
library(readr)
library(tidyverse)
library(stats)

# Read the data
data <- read_csv("data_files/sample_points_in_10D.csv")

# Set row names and remove the first column
rownames(data) <- data$...1
data <- data[ , -1]

# Calculate pairwise distances in the original 10-dimensional space
distances_original <- dist(data)



# Perform PCA to reduce dimensions to 2
pca_result <- prcomp(data, scale. = TRUE)
data_pca <- as.data.frame(pca_result$x[, 1:2])

# Calculate pairwise distances in the reduced 2-dimensional space
distances_reduced <- dist(data_pca)

# Print the pairwise distances
print("Pairwise distances in the original 10-dimensional space:")
print(as.matrix(distances_original))
hist(distances_original)

print("Pairwise distances in the reduced 2-dimensional space:")
print(as.matrix(distances_reduced))
hist(distances_reduced)


# Convert distance objects to data frames
distances_original_df <- as.data.frame(as.matrix(distances_original))
distances_reduced_df <- as.data.frame(as.matrix(distances_reduced))

# Convert data frames to long format for ggplot2
distances_original_long <- distances_original_df %>%
  rownames_to_column("Point1") %>%
  gather(Point2, Distance_Original, -Point1)

distances_reduced_long <- distances_reduced_df %>%
  rownames_to_column("Point1") %>%
  gather(Point2, Distance_Reduced, -Point1)

# Merge the two data frames on Point1 and Point2
distances_comparison <- merge(distances_original_long, distances_reduced_long, by = c("Point1", "Point2"))


p <- ggplot(distances_comparison, aes(x = Distance_Original, y = Distance_Reduced)) +
  geom_point(color = "blue", size = 0.95) +
  geom_abline(slope = 1, intercept = 0, color = "green", linetype = "dashed") +
  annotate("text", 
           x = max(distances_comparison$Distance_Original) * 0.2,
           y = max(distances_comparison$Distance_Original) * 0.2,
           label = "y = x",
           angle = 17,
           vjust = -0.5) +
  labs(x = "Distances in 10D space",
       y = "Distances in 2D space") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 8),  # Set axis labels to size 8
    axis.text = element_text(size = 8)    # Set tick numbers to size 8
  )

ggsave_white_bg("pairwise_distance_10d_2d_comparison.pdf", 
       plot = p,
       width = 3.5, 
       height = 3.5,
       units = "in",
       dpi = 300)
