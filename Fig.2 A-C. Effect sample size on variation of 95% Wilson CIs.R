# Load necessary packages
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(tidyr)


# Parameters
set.seed(500910000)  # For reproducibility
Prevalence_vals <- c(5, 15, 25, 35, 50, 65, 75, 85, 95) # Given prevalence values
B <- 1000            # Number of little bootstrap samples
b <- 100             # Bootstrap resamples for each little sample
n_vals <- seq(5, 300, by = 5) # Subsample sizes

# Perform Bag of Little Bootstraps (BLB) for MEDIAN
perform_blb <- function(Prevalence, n_vals, B, b) {
  results <- list()
  for (n in n_vals) {
    estimates <- numeric(B)
    for (i in 1:B) {
      # 1. Create a bootstrap sample of size n
      little_sample <- rbinom(n, 1, Prevalence / 100)
      # 2. Apply the ordinary bootstrap to the bootstrap sample
      boot_samples <- numeric(b)
      for (j in 1:b) {
        resample_indices_boot <- sample(1:n, n, replace = TRUE)
        resample_boot <- little_sample[resample_indices_boot]
        boot_samples[j] <- mean(resample_boot) * 100 # Convert to percentage
      }
      # 3. Calculate the median prevalence of the bootstrap samples
      estimates[i] <- median(boot_samples)
    }
    results[[as.character(n)]] <- estimates
  }
  return(results)
}

# Calculate Wilson score interval and precision
calculate_wilson_ci_precision <- function(Prevalence, n_vals) {
  results <- list()
  for (n in n_vals) {
    # Calculate Wilson score interval
    p <- Prevalence / 100
    z <- 1.96  # Z-score for 95% CI
    lower_wilson <- (p + z^2 / (2 * n) - z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n)) / (1 + z^2 / n) * 100
    upper_wilson <- (p + z^2 / (2 * n) + z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n)) / (1 + z^2 / n) * 100
    precision <- (upper_wilson - lower_wilson) / 2 # Calculate precision
    results[[as.character(n)]] <- list(ci_lower = lower_wilson, ci_upper = upper_wilson, precision = precision)
  }
  return(results)
}

# Apply BLB to each prevalence value
all_blb_results <- list()
for (Prevalence in Prevalence_vals) {
  all_blb_results[[as.character(Prevalence)]] <- perform_blb(Prevalence, n_vals, B, b)
}

# Calculate Wilson CI and precision for each prevalence
all_wilson_results <- list()
for (Prevalence in Prevalence_vals) {
  all_wilson_results[[as.character(Prevalence)]] <- calculate_wilson_ci_precision(Prevalence, n_vals)
}

# --- Data Preparation for Plotting ---
plot_data <- data.frame(
  n = numeric(),
  Prevalence = numeric(),
  Median = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  Precision = numeric(),
  Prevalence = character(),
  stringsAsFactors = FALSE
)

for (Prevalence in Prevalence_vals) {
  for (n in n_vals) {
    n_char <- as.character(n)
    # Extract values
    median_val <- median(all_blb_results[[as.character(Prevalence)]][[n_char]], na.rm = TRUE)
    ci_lower_val <- all_wilson_results[[as.character(Prevalence)]][[n_char]]$ci_lower
    ci_upper_val <- all_wilson_results[[as.character(Prevalence)]][[n_char]]$ci_upper
    precision_val <- all_wilson_results[[as.character(Prevalence)]][[n_char]]$precision
    temp_data <- data.frame(
      n = n,
      Prevalence = Prevalence,
      Median = median_val,
      CI_Lower = ci_lower_val,
      CI_Upper = ci_upper_val,
      Precision = precision_val,
      Expected_prevalence = paste0(Prevalence, "%")
    )
    plot_data <- rbind(plot_data, temp_data)
  }
}

plot_data <- plot_data %>% filter(!is.na(Median))
plot_data$Expected_prevalence <- factor(plot_data$Expected_prevalence, levels = unique(plot_data$Expected_prevalence))
color_palette <- brewer.pal(n = 9, name = "Set1")

# Plot 1: Median prevalence and confidence intervals
plot_1 <-ggplot(plot_data, aes(x = n, y = Median, color = Expected_prevalence, group = Expected_prevalence)) +
  geom_hline(aes(yintercept = Prevalence, color = Expected_prevalence),
             linetype = "solid", linewidth = 0.5) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper, fill = Expected_prevalence, color = Expected_prevalence),
              alpha = 0.25,
              linetype = "dotted",
              linewidth = 0.2) +
  # Add vertical lines for n=25, 100 and 150
  annotate("segment", x = 15, xend = 15, y = 20, yend = 80,
           linewidth = 0.7, color = "black")  +
   annotate("segment", x = 85, xend = 85, y = 2, yend = 10,
           linewidth = 0.7, color = "black")+  
  annotate("segment", x = 85, xend = 85, y = 90, yend = 98,
           linewidth = 0.7, color = "black")+ 
 annotate("segment", x = 25, xend = 25, y = 10, yend = 20,
           linewidth = 0.7, color = "black")+
  annotate("segment", x = 25, xend = 25, y = 80, yend = 90,
           linewidth = 0.7, color = "black")+
  geom_vline(xintercept = 110, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 135, color = "blue", linetype = "dashed") +
  annotate("rect", xmin = 111, xmax = 134, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.25) +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  scale_x_continuous(breaks = c(5, 20, 40, 60, 80, 100, 125, 150, 175, 200),
                     labels = as.character(c(5, 20, 40, 60, 80, 100, 125, 150, 175, 200)),
                     limits = c(0, 200)) +
  scale_y_continuous(breaks = c(5, 15,  25,  35,  45,  55,  65,  75,  85,  95),
                     labels = as.character(c(5, 15,  25,  35,  45,  55,  65,  75,  85,  95)),
                     limits = c(0, 100))+
  geom_line(linetype = "dashed") +
  geom_point(size = 1) +
  labs(
    x = "Sample size",
    y = "Prevalence (%)",
    tag = "A" # Add tag
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Plot 2: Precision
# plot_data <- data.frame(...) 
# Filter the data for the specified prevalences
filtered_data <- plot_data %>%
  filter(Expected_prevalence %in% c("5%", "15%", "25%", "35%", "50%"))

# Create a named vector for the new labels
new_labels <- c(
  "5%" = "5% / 95%",
  "15%" = "15% / 85%",
  "25%" = "25% / 75%",
  "35%" = "35% / 65%",
  "50%" = "50%"
)

# Apply the new labels to the Expected_prevalence column
filtered_data <- filtered_data %>%
  mutate(Expected_prevalence = recode(Expected_prevalence, !!!new_labels))

# Ensure Expected_prevalence is a factor with the new labels as levels
filtered_data$Expected_prevalence <- factor(filtered_data$Expected_prevalence, levels = new_labels[unique(filtered_data$Expected_prevalence)])

# Now filtered_data contains the filtered data with the new labels

plot_2 <-ggplot(filtered_data,aes(x = n, y = Precision, color = Expected_prevalence, group = Expected_prevalence)) +
  geom_line() +
  geom_point(size = 0.9) +
  scale_color_manual(values = c(
    "5% / 95%" = "red",       # Example color for 5% / 95%
    "15% / 85%" = "blue",      # Example color for 15% / 85
    "25% / 75%" = "green",     # Example color for 25% / 75%
    "35% / 65%" = "yellow3",    # Example color for 35% / 65%
    "50%"      = "brown4"     # Example color for 50%
  )) +
  scale_x_continuous(breaks = c(5, 20, 40, 60, 80, 100, 125, 150, 175, 200, 225, 250, 275, 300),
                     labels = as.character(c(5, 20, 40, 60, 80, 100, 125, 150, 175, 200, 225, 250, 275, 300)),
                     limits = c(5, 300)) +
  scale_y_continuous(breaks = c(2, 5, 10, 15, 20, 25, 30, 35),
                     labels = as.character(c(2, 5, 10, 15, 20, 25, 30, 35)),
                     limits = c(2, 35)) +
  labs(
    x = "Sample size",
    y = "Precision (%)",
    color = "Expected_prevalence",  # Changed the label to "Prevalence"
    tag = "B"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))




# Dynamically generate ranges
increment <- n_vals[2] - n_vals[1]
ranges <- list()
for (i in 1:(length(n_vals) - 1)) {
  start <- n_vals[i]
  end <- n_vals[i+1]
  range_name <- paste0(start, "-", end)
  ranges[[range_name]] <- start:end
}


calculate_precision_change <- function(data) { # Removed range_values argument as it's not used
  if (nrow(data) < 2) return(NA)
  precision_change <- max(data$Precision, na.rm = TRUE) - min(data$Precision, na.rm = TRUE) # Added na.rm = TRUE
  return(precision_change)
}

precision_changes <- do.call(rbind, lapply(unique(plot_data$Expected_prevalence), function(region) {
  region_data <- plot_data %>% filter(Expected_prevalence == region)

  # Apply calculate_precision_change to each range and handle potential errors
  precision_change_values <- sapply(ranges, function(range_values) {
    subset_data <- filter(region_data, n %in% range_values)
    if (nrow(subset_data) == 0) {
      return(NA) # Return NA if subset is empty
    }
    calculate_precision_change(subset_data) # Pass the subsetted data
  })

  data.frame(
    Expected_prevalence = region,
    Range = names(ranges),
    PrecisionChange = precision_change_values # Use the calculated values
  )
}))

precision_changes$Range <- factor(precision_changes$Range, levels = names(ranges))

#Corrected object name here and removed the %>% operator, as it's not needed
filtered_data_ch <- precision_changes %>%
  filter(Expected_prevalence %in% c("5%", "15%", "25%", "35%", "50%"))

# Create a named vector for the new labels
new_labels <- c(
 "5%" = "5% / 95%",
  "15%" = "15% / 85%",
  "25%" = "25% / 75%",
  "35%" = "35% / 65%",
  "50%" = "50%"
)

# Apply the new labels to the Expected_prevalence column
filtered_data_ch <- filtered_data_ch %>%  # Corrected object name here
  mutate(Expected_prevalence = recode(Expected_prevalence, !!!new_labels))

# Ensure Expected_prevalence is a factor with the new labels as levels
filtered_data_ch$Expected_prevalence <- factor(filtered_data_ch$Expected_prevalence, levels = new_labels[unique(filtered_data_ch$Expected_prevalence)])


plot_3 <-ggplot(filtered_data_ch, aes(x = Range, y = PrecisionChange, color = Expected_prevalence, group = Expected_prevalence)) +
  geom_line() + 
  geom_point(size = 1) +
  scale_color_manual(values = c(
    "5% / 95%" = "red",       # Example color for 5% / 95%
    "15% / 85%" = "blue",      # Example color for 15% / 85
    "25% / 75%" = "green",     # Example color for 25% / 75%
    "35% / 65%" = "yellow3",    # Example color for 35% / 65%
    "50%"      = "brown4"     # Example color for 50%
  )) +
  labs(x = "Sample size increment", y = "Change in precision", 
  Color = "Expected_prevalence", tag = "C") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "false")+
  scale_x_discrete(breaks = levels(precision_changes$Range)[seq(1, length(levels(precision_changes$Range)), by = 5)])
  

# Combine plots
combined_plot <- plot_1 + plot_2 + plot_3 + plot_layout(guides = "collect")
print(combined_plot)

#Export data in CSV
write.csv(filtered_data, file = "Data to fig 2A.csv", row.names = T)


