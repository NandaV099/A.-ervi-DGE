
########################### Housekeeping ##############################
# Set working directory
setwd("C:/Users/nanda/Desktop/Master Thesis/Morphometrics/TPS")

rm(list=ls()) 




# Load Libraries
library(ggfortify)
library(FactoMineR)
library(lattice)
library(vcd)
library(factoextra)
library(stats)
library(plyr)
library(dplyr)
library(ggcorrplot)
library(corrr)
library(zoom)            
library(geomorph)
library(abind) #combines multidimensional data frames like my tps files
library(ggplot2)
library(ggpubr)
library(plotrix)
library(vegan)
library(Morpho)
library(finalfit)
library(MASS)
library(rstan)
library(data.table)
library(sjPlot)
library(boot)
library(stargazer)
library(caret) 
library(lme4)
library(FactoMineR)
library(ggfortify)



# Load required libraries
library(dplyr)




########################### Housekeeping ##############################
# Set working directory
setwd("C:/Users/nanda/Desktop/Master Thesis/Morphometrics/TPS")

rm(list = ls())  # Clear the workspace




########################## Load and Parse TPS Files #################################
# List of TPS files for Wings and Legs
tps_files_wings <- list(
  "FA_FA_Wing_F34.TPS",
  "FA_PA_Wing_F34.TPS",
  "PA_FA_Wing_F34.TPS",
  "PA_PA_Wing_F34.TPS",
  "FA_FA_Wing_F39.TPS",
  "FA_PA_Wing_F39.TPS",
  "PA_FA_Wing_F39.TPS",
  "PA_PA_Wing_F39.TPS"
)

tps_files_legs <- list(
  "FA_FA_Leg_F34.TPS",
  "FA_PA_Leg_F34.TPS",
  "PA_FA_Leg_F34.TPS",
  "PA_PA_Leg_F34.TPS",
  "FA_FA_Leg_F39.TPS",
  "FA_PA_Leg_F39.TPS",
  "PA_FA_Leg_F39.TPS",
  "PA_PA_Leg_F39.TPS"
)

# Function to manually parse TPS files and retain correct IDs and coordinates
parse_tps_file <- function(tps_file_path) {
  lines <- readLines(tps_file_path)
  landmarks <- list()
  current_landmarks <- NULL
  current_id <- NULL
  scale_factor <- NA
  
  # Determine the family based on the file name
  family_tag <- ifelse(grepl("F34", tps_file_path), "F34", "F39")
  
  for (line in lines) {
    if (grepl("^LM=", line)) {
      # Start new set of landmarks
      if (!is.null(current_landmarks)) {
        current_id <- sub("_F\\d{2}_", paste0("_", family_tag, "_"), current_id)
        landmarks[[current_id]] <- list(landmarks = current_landmarks, scale = scale_factor)
      }
      current_landmarks <- matrix(ncol = 2, nrow = as.numeric(sub("LM=", "", line)))
    } else if (grepl("^\\d+\\.\\d+", line)) {
      # Extract coordinates
      coords <- as.numeric(unlist(strsplit(line, "\\s+")))
      current_landmarks[which(is.na(current_landmarks[, 1]))[1], ] <- coords
    } else if (grepl("^ID=", line)) {
      # Set specimen ID without modifying it
      current_id <- sub("ID=", "", line)
    } else if (grepl("^SCALE=", line)) {
      # Set scale factor
      scale_factor <- as.numeric(sub("SCALE=", "", line))
    }
  }
  
  # Add the last set of landmarks
  if (!is.null(current_landmarks) && !is.null(current_id)) {
    current_id <- sub("_F\\d{2}_", paste0("_", family_tag, "_"), current_id)
    landmarks[[current_id]] <- list(landmarks = current_landmarks, scale = scale_factor)
  }
  
  return(landmarks)
}

# Parse all TPS files for Wings
landmarks_data_wings <- list()
for (file in tps_files_wings) {
  file_path <- paste0("C:/Users/nanda/Desktop/Master Thesis/Morphometrics/TPS/", file)
  parsed_data <- parse_tps_file(file_path)
  landmarks_data_wings <- append(landmarks_data_wings, parsed_data)
}

# Parse all TPS files for Legs
landmarks_data_legs <- list()
for (file in tps_files_legs) {
  file_path <- paste0("C:/Users/nanda/Desktop/Master Thesis/Morphometrics/TPS/", file)
  parsed_data <- parse_tps_file(file_path)
  landmarks_data_legs <- append(landmarks_data_legs, parsed_data)
}




########################## Create Data Frames for Wings and Legs #############################

# Extract the dimension names from wings and legs
dim_names_wings <- names(landmarks_data_wings)
dim_names_legs <- names(landmarks_data_legs)

# Initialize empty vectors to store metadata for wings and legs
create_metadata <- function(dim_names) {
  Line <- vector()
  Conditioned <- vector()
  Family <- vector()
  
  # Iterate through each name in dim_names to extract Line, Conditioned, and Family
  for (name in dim_names) {
    # Split the name by underscore
    parts <- unlist(strsplit(name, "_"))
    
    # Extract information based on the naming convention
    Line <- c(Line, parts[1]) # The first part (FA or PA)
    Conditioned <- c(Conditioned, parts[2]) # The second part (FA or PA)
    Family <- c(Family, parts[3]) # The family (F34 or F39)
  }
  
  # Create a data frame to hold metadata
  metadata <- data.frame(Line = Line, Conditioned = Conditioned, Family = Family)
  return(metadata)
}

metadata_wings <- create_metadata(dim_names_wings)
metadata_legs <- create_metadata(dim_names_legs)




########################## Calculate Wing Surface Area #############################

# Function to calculate the surface area of wings using the shoelace formula
calculate_surface_area_manual <- function(landmark_matrix, scale_factor) {
  # Extract x and y coordinates for each landmark
  x_coords <- landmark_matrix[, 1]
  y_coords <- landmark_matrix[, 2]
  
  # Use the shoelace formula to calculate the area of the polygon
  area <- 0.5 * abs(sum(x_coords * c(y_coords[-1], y_coords[1]) - y_coords * c(x_coords[-1], x_coords[1])))
  
  # Apply scale factor to convert the area to micrometers squared
  if (is.na(scale_factor)) {
    warning("Scale factor is NA, returning area without scaling.")
    scaled_area <- area
  } else {
    scaled_area <- area * (scale_factor ^ 2)
  }
  
  return(scaled_area)
}

# Calculate Wing Surface Areas
wing_surface_areas <- numeric(length(landmarks_data_wings))
for (i in 1:length(landmarks_data_wings)) {
  landmark_matrix <- landmarks_data_wings[[i]]$landmarks
  scale_factor <- landmarks_data_wings[[i]]$scale
  wing_surface_areas[i] <- calculate_surface_area_manual(landmark_matrix, scale_factor)
}

# Add Wing Surface Area to Metadata
metadata_wings$Wing_Surface_Area_um2 <- wing_surface_areas




########################## Calculate Tibia Length #############################

# Function to calculate tibia length as distance between landmark 1 and landmark 2
calculate_tibia_length <- function(landmark_matrix, scale_factor) {
  # Extract x and y coordinates for landmark 1 and 2
  x1 <- landmark_matrix[1, 1]
  y1 <- landmark_matrix[1, 2]
  x2 <- landmark_matrix[2, 1]
  y2 <- landmark_matrix[2, 2]
  
  # Calculate Euclidean distance
  length <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  
  # Apply scale factor to convert to micrometers
  if (is.na(scale_factor)) {
    warning("Scale factor is NA, returning length without scaling.")
    scaled_length <- length
  } else {
    scaled_length <- length * scale_factor
  }
  
  return(scaled_length)
}

# Calculate Tibia Lengths
tibia_lengths <- numeric(length(landmarks_data_legs))
for (i in 1:length(landmarks_data_legs)) {
  landmark_matrix <- landmarks_data_legs[[i]]$landmarks
  scale_factor <- landmarks_data_legs[[i]]$scale
  tibia_lengths[i] <- calculate_tibia_length(landmark_matrix, scale_factor)
}

# Add Tibia Length to Metadata
metadata_legs$Tibia_Length_um <- tibia_lengths




######################### Tibia Length Analysis (Legs) #########################

# Assuming `metadata_legs` contains the tibia length measurements

# Linear models for Tibia Length
model_leg_1 <- lm(`Tibia_Length_um` ~ Line + Conditioned + Family, data = metadata_legs)
model_leg_2 <- lm(`Tibia_Length_um` ~ Line * Conditioned + Family, data = metadata_legs)
model_leg_3 <- lm(`Tibia_Length_um` ~ Line * Family + Conditioned, data = metadata_legs)
model_leg_4 <- lm(`Tibia_Length_um` ~ Line + Conditioned * Family, data = metadata_legs)
model_leg_5 <- lmer(Tibia_Length_um ~ Line * Conditioned + (1 | Family), data = metadata_legs)

# Summary of models for Tibia Length
summary(model_leg_1)
summary(model_leg_2)
summary(model_leg_3)
summary(model_leg_4)
summary(model_leg_5)


# AIC comparison between leg models
aic_values_leg <- AIC(model_leg_1, model_leg_2, model_leg_3, model_leg_4, model_leg_5)
print(aic_values_leg)

# Prepare data for plotting Tibia Length
metadata_legs$Treatment <- paste(metadata_legs$Line, metadata_legs$Conditioned, sep = "_")

# Filter data for Family F34 and F39 (Legs)
leg_f34_data <- subset(metadata_legs, Family == "F34")
leg_f39_data <- subset(metadata_legs, Family == "F39")

# Boxplots for Tibia Length for F34 and F39
plot_leg_f34 <- ggplot(leg_f34_data, aes(x = Treatment, y = `Tibia_Length_um`, fill = Treatment)) +
  geom_boxplot() +
  labs(title = "Tibia Length Comparison for F34 (Legs)", x = "Line - Conditioned", y = "Tibia Length (µm)") +
  theme_minimal() +
  theme(legend.position = "none")

plot_leg_f39 <- ggplot(leg_f39_data, aes(x = Treatment, y = `Tibia_Length_um`, fill = Treatment)) +
  geom_boxplot() +
  labs(title = "Tibia Length Comparison for F39 (Legs)", x = "Line - Conditioned", y = "Tibia Length (µm)") +
  theme_minimal() +
  theme(legend.position = "none")

# Arrange and display the plots side by side for Tibia Length
combined_plot_leg <- ggarrange(plot_leg_f34, plot_leg_f39, ncol = 2, nrow = 1)
print(combined_plot_leg)

# ANOVA and Tukey's HSD for Legs
aov_leg_f34 <- aov(`Tibia_Length_um` ~ Line * Conditioned, data = leg_f34_data)
summary(aov_leg_f34)
tukey_leg_f34 <- TukeyHSD(aov_leg_f34)
print(tukey_leg_f34)

aov_leg_f39 <- aov(`Tibia_Length_um` ~ Line * Conditioned, data = leg_f39_data)
summary(aov_leg_f39)
tukey_leg_f39 <- TukeyHSD(aov_leg_f39)
print(tukey_leg_f39)

# T-tests between F34 and F39 for each treatment (Legs)
treatments_leg <- unique(metadata_legs$Treatment)
for (treatment in treatments_leg) {
  leg_f34_subset <- subset(metadata_legs, Treatment == treatment & Family == "F34")
  leg_f39_subset <- subset(metadata_legs, Treatment == treatment & Family == "F39")
  
  # Perform t-test if there are enough samples
  if (nrow(leg_f34_subset) > 1 && nrow(leg_f39_subset) > 1) {
    t_test_leg_result <- t.test(leg_f34_subset$`Tibia_Length_um`, leg_f39_subset$`Tibia_Length_um`)
    cat(sprintf("\nT-test for Treatment %s (F34 vs F39, Legs):\n", treatment))
    print(t_test_leg_result)
  } else {
    cat(sprintf("\nNot enough data for t-test for Treatment %s (F34 vs F39, Legs).\n", treatment))
  }
}


######################### Wing Surface Area Analysis (Wings) #########################
# Load necessary libraries
library(ggplot2)
library(viridis)

# Add unique specimen identifiers to each dataframe for merging
metadata_legs$Specimen_ID <- seq_len(nrow(metadata_legs))
metadata_wings$Specimen_ID <- seq_len(nrow(metadata_wings))

# Merge datasets based on the unique identifier
merged_data <- merge(metadata_legs, metadata_wings, by = "Specimen_ID", suffixes = c("_legs", "_wings"))

# Create a new column for Sub-Treatment (Line - Conditioned combination)
merged_data$Sub_Treatment <- interaction(merged_data$Line_legs, merged_data$Conditioned_legs, sep = " - ")

# Journal of Evolutionary Biology theme (JEB-friendly)
jeb_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_blank(), # Remove box border
    panel.grid.major = element_blank(), # Remove gridlines
    panel.grid.minor = element_blank(), # Remove minor gridlines
    axis.line = element_line(size = 0.8, color = "black"), # X and Y lines only
    legend.position = "bottom", # Legend at the bottom
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  )

# Colorblind-friendly palette
colorblind_palette <- c(
  "FA - FA" = "#0072B2", # Blue
  "FA - PA" = "#E69F00", # Orange
  "PA - FA" = "#009E73", # Green
  "PA - PA" = "#CC79A7"  # Pink
)

# Calculate correlation coefficient and p-value
cor_test_result <- cor.test(merged_data$Tibia_Length_um, merged_data$Wing_Surface_Area_mm2)

# Extract values for display in the annotation
cor_coeff <- round(cor_test_result$estimate, 2)
p_value <- signif(cor_test_result$p.value, 2)

# Plot with JEB formatting and colorblind-friendly palette
correlation_plot <- ggplot(merged_data, aes(x = Tibia_Length_um, y = Wing_Surface_Area_mm2, color = Sub_Treatment)) +
  geom_point(size = 2, alpha = 0.8) +  # Grayscale-friendly points
  geom_smooth(method = "lm", color = "black", size = 1, se = TRUE) +  # Black regression line
  scale_color_manual(values = colorblind_palette, name = "Line - Conditioned:") +
  labs(
    title = "Correlation between Tibia Length and Wing Surface Area",
    x = "Tibia Length (µm)",
    y = "Wing Surface Area (mm²)"
  ) +
  annotate(
    "text", x = min(merged_data$Tibia_Length_um, na.rm = TRUE), 
    y = max(merged_data$Wing_Surface_Area_mm2, na.rm = TRUE) - 0.05 * diff(range(merged_data$Wing_Surface_Area_mm2, na.rm = TRUE)), 
    label = paste("r =", cor_coeff, "\np =", p_value),
    hjust = 0, vjust = 1, size = 5, fontface = "bold"
  ) +
  jeb_theme

# Save the plot
ggsave(
  "Correlation_Plot_Colorblind_Friendly.png",
  plot = correlation_plot,
  bg = "white",  # Ensure white background
  width = 8,     # JEB-suitable width
  height = 6,    # JEB-suitable height
  dpi = 300
)

# Display the plot
print(correlation_plot)


######################### Combined plots #########################
# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(viridis)

# Ensure `Treatment` column matches the scale fill values
metadata_legs$Treatment <- gsub("_", " - ", metadata_legs$Treatment)
metadata_wings$Treatment <- gsub("_", " - ", metadata_wings$Treatment)

# Define the treatment colors using the viridis palette (colorblind-friendly)
treatment_colors <- c(
  "FA - FA" = "#E69F00",  # Orange
  "FA - PA" = "#56B4E9",  # Blue
  "PA - FA" = "#009E73",  # Green
  "PA - PA" = "#CC79A7"   # Purple
)

# Check the unique levels in Treatment
print("Unique levels in metadata_legs:")
print(unique(metadata_legs$Treatment))
print("Unique levels in metadata_wings:")
print(unique(metadata_wings$Treatment))

# JEB-style theme (based on the example you provided)
jeb_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Smaller font for main header
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    axis.line = element_line(size = 0.8, color = "black"),  # Thicker axis lines
    axis.ticks = element_line(size = 0.8, color = "black"),  # Axis ticks
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),  # No major gridlines
    panel.grid.minor = element_blank(),  # No minor gridlines
    panel.border = element_blank(),  # No panel border
    axis.text.x = element_blank(),  # No x-axis labels
    axis.ticks.x = element_blank()  # No x-axis ticks
  )

# Create the first plot (Tibia Length Comparison)
plot_leg_combined <- ggplot(metadata_legs, aes(x = Treatment, y = `Tibia_Length_um`, fill = Treatment)) +
  geom_boxplot(color = "black", linewidth = 0.8) +
  facet_wrap(~ Family) +
  labs(
    x = NULL, 
    y = "Tibia Length (µm)",
    title = "A) Tibia Length Comparison"
  ) +
  scale_fill_viridis_d(option = "C", direction = 1, name = "Line - Conditioned") +  # Viridis color palette
  jeb_theme

# Create the second plot (Wing Surface Area Comparison)
plot_wing_combined <- ggplot(metadata_wings, aes(x = Treatment, y = Wing_Surface_Area_mm2, fill = Treatment)) +
  geom_boxplot(color = "black", linewidth = 0.8) +
  facet_wrap(~ Family) +
  labs(
    x = NULL, 
    y = "Wing Surface Area (mm²)",
    title = "B) Wing Surface Area Comparison"
  ) +
  scale_fill_viridis_d(option = "C", direction = 1, name = "Line - Conditioned") +  # Viridis color palette
  jeb_theme

# Combine the plots with a shared legend
final_combined_plot <- ggarrange(
  plot_leg_combined, 
  plot_wing_combined, 
  ncol = 2, 
  common.legend = TRUE, 
  legend = "bottom"
)

# Display the plot
print(final_combined_plot)

# Save the final combined plot
ggsave(
  "Tibia+Wing_Final.png",
  final_combined_plot,
  bg = "white",  # Force white background
  width = 12,
  height = 6,
  units = "in",
  dpi = 300
)




######################### PCA for tibia length and wing surface area #########################
# Load required libraries
library(ggplot2)
library(viridis)
library(gridExtra)
library(geomorph)

# Define color palette for 'Line - Conditioned' combinations
treatment_colors <- c(
  "FA - FA" = "#E69F00",  # Orange
  "FA - PA" = "#56B4E9",  # Blue
  "PA - FA" = "#009E73",  # Green
  "PA - PA" = "#CC79A7"   # Purple
)

# Step 1: Perform Generalized Procrustes Analysis (GPA) if not already done
wing_coords_array <- array(
  unlist(lapply(landmarks_data_wings, function(x) x$landmarks)),
  dim = c(nrow(landmarks_data_wings[[1]]$landmarks), 2, length(landmarks_data_wings))
)
gpa_wings <- gpagen(wing_coords_array)

# Step 2: Combine the aligned coordinates with metadata in a geomorph data frame
all_info_wings <- geomorph.data.frame(
  coords = gpa_wings$coords,
  Line = metadata_wings$Line,
  Conditioned = metadata_wings$Conditioned,
  Family = metadata_wings$Family,
  Centroid_Size = gpa_wings$Csize,
  Wing_Surface_Area_mm2 = metadata_wings$Wing_Surface_Area_mm2
)

# Step 3: Run PCA on the aligned coordinates
pc_wings <- gm.prcomp(all_info_wings$coords)

# Step 4: Create a data frame for visualization with centroid size and surface area included
pca_data_wings <- data.frame(
  pc1 = pc_wings$x[, 1],
  pc2 = pc_wings$x[, 2],
  pc3 = pc_wings$x[, 3],
  Centroid_Size = all_info_wings$Centroid_Size,
  Wing_Surface_Area_mm2 = all_info_wings$Wing_Surface_Area_mm2,
  Line = all_info_wings$Line,
  Conditioned = all_info_wings$Conditioned,
  Family = all_info_wings$Family
)

# Add a combined treatment variable for plotting
pca_data_wings$Treatment <- interaction(pca_data_wings$Line, pca_data_wings$Conditioned, sep = " - ")

# Step 5: Combine centroid size and surface area for PCA
pca_data_with_size_proxies <- data.frame(
  Centroid_Size = pca_data_wings$Centroid_Size,
  Wing_Surface_Area_mm2 = pca_data_wings$Wing_Surface_Area_mm2,
  pc1 = pca_data_wings$pc1,
  pc2 = pca_data_wings$pc2
)

# Perform PCA including both centroid size and surface area
pca_combined_size <- prcomp(pca_data_with_size_proxies, scale. = TRUE)

# Create a data frame for plotting the PCA results with both size proxies included
pca_results_combined_size <- data.frame(
  PC1 = pca_combined_size$x[, 1],
  PC2 = pca_combined_size$x[, 2],
  Treatment = pca_data_wings$Treatment,
  Line = pca_data_wings$Line,
  Conditioned = pca_data_wings$Conditioned,
  Family = pca_data_wings$Family
)

# Update the 'Treatment' variable to ensure consistent color mapping
pca_results_combined_size$Treatment <- factor(
  pca_results_combined_size$Treatment, 
  levels = c("PA - PA", "FA - PA", "PA - FA", "FA - FA")
)

# Calculate the explained variance for PC1 and PC2
explained_variance <- summary(pca_combined_size)$importance[2, ]
pc1_var <- round(explained_variance[1] * 100, 2)
pc2_var <- round(explained_variance[2] * 100, 2)

# Plot PCA by Line - Conditioned combinations
pca_plot_conditioned <- ggplot(pca_results_combined_size, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 3, alpha = 1) +  # Fully opaque points
  stat_ellipse(aes(color = Treatment), alpha = 1, linetype = "solid", size = 1.2, level = 0.95) +  # Ellipses
  scale_color_manual(values = treatment_colors) +
  xlab(paste0("PC1 (", pc1_var, "% Variance)")) +
  ylab(paste0("PC2 (", pc2_var, "% Variance)")) +
  ggtitle("Wing Shape by Treatment") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.key = element_blank(),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(size = 0.8, color = "black")
  )

# Plot PCA by Family (ignoring Treatment)
pca_plot_family <- ggplot(pca_results_combined_size, aes(x = PC1, y = PC2, color = Family)) +
  geom_point(size = 3, alpha = 1) +  # Fully opaque points
  stat_ellipse(aes(color = Family), alpha = 1, linetype = "solid", size = 1.2, level = 0.95) +  # Ellipses
  scale_color_manual(values = c("F34" = "#F4A300", "F39" = "#1E78B5")) +  # Custom colors for families
  xlab(paste0("PC1 (", pc1_var, "% Variance)")) +
  ylab(paste0("PC2 (", pc2_var, "% Variance)")) +
  ggtitle("Wing Shape by Family") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.key = element_blank(),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(size = 0.8, color = "black")
  )

# Combine both PCA plots into a single figure
final_combined_plot <- ggarrange(
  pca_plot_conditioned,
  pca_plot_family,
  ncol = 2, nrow = 1,
  common.legend = FALSE
)

# Display the plot
print(final_combined_plot)

# Save the combined plot
ggsave(
  "PCA_Plot_Final.png",
  final_combined_plot,
  width = 12,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
) 
