# Load necessary packages
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(tidyr)

# Empirical data on the prevalence of Echinococcus multilocularis in red foxes from eight regions in Slovakia were accumulated over 5 years (2000, 2001, 2002, 2005/2006, and 2007), as reported by Miterpáková and Dubinský (2011) DOI: 10.2478/s11687-011-0023-5
empirical_data <- data.frame(
  Region = c("BA", "TT", "NR", "TN", "ZA", "BB", "KE", "PO"),
  dissected = c(172, 271, 476, 254, 379, 762, 670, 647),
  infected = c(22, 42, 118, 105, 197, 228, 127, 259),
  proportion = c(0.13, 0.15, 0.25, 0.41, 0.52, 0.30, 0.19, 0.40)
)

# Parameters
set.seed(500910000)
B <- 1000  # Number of bootstrap samples
b <- 100   # Bootstrap resamples within each bootstrap sample
n_vals <- c(3, 6, 9, 12, 14, 19, 24, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300, 305, 310, 315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390, 395, 400, 405, 410, 415, 420, 425, 430, 435, 440, 445, 450) # Subsample sizes

# Perform Bag of Little Bootstraps (BLB) - for MEDIAN
perform_blb_empirical <- function(data_row, n_vals, B, b) {
  results <- list()

  for (n in n_vals) {
    estimates <- numeric(B)

    for (i in 1:B) {
      # 1. Create a bootstrap sample of size n
      if (n > data_row$dissected) {
        estimates[i] <- NA  # Return NA if n exceeds dissected count
        next
      }

      resample_indices <- sample(1:data_row$dissected, n, replace = TRUE)

      # Create a bootstrap sample based on the resampled indices
      little_sample <- ifelse(resample_indices <= data_row$infected, 1, 0)

      # 2. Apply the ordinary bootstrap to the bootstrap sample
      boot_samples <- numeric(b)
      for (j in 1:b) {
        # Resample with replacement from the bootstrap sample
        resample_indices_boot <- sample(1:n, n, replace = TRUE)
        resample_boot <- little_sample[resample_indices_boot]

        # Calculate prevalence in the resample
        boot_samples[j] <- mean(resample_boot) * 100
      }

      # 3. Calculate the median prevalence of the bootstrap samples
      estimates[i] <- median(boot_samples)
    }
    results[[as.character(n)]] <- estimates
  }
  return(results)
}

# Calculate Wilson score interval and precision - from EMPIRICAL DATA
calculate_wilson_ci_precision <- function(data_row, n_vals) {
  results <- list()

  for (n in n_vals) {
    if (n > data_row$dissected) {
      results[[as.character(n)]] <- list(ci_lower = NA, ci_upper = NA, precision = NA)
      next
    }

    # Calculate Wilson score interval
    p <- data_row$infected / data_row$dissected
    z <- 1.96  # Z-score for 95% CI
    n_effective <- n  # Use the actual sample size

    lower_wilson <- (p + z^2 / (2 * n_effective) - z * sqrt((p * (1 - p) + z^2 / (4 * n_effective)) / n_effective)) / (1 + z^2 / n_effective) * 100
    upper_wilson <- (p + z^2 / (2 * n_effective) + z * sqrt((p * (1 - p) + z^2 / (4 * n_effective)) / n_effective)) / (1 + z^2 / n_effective) * 100

    # Calculate precision (width of the CI)
    precision <- (upper_wilson - lower_wilson) / 2

    results[[as.character(n)]] <- list(ci_lower = lower_wilson, ci_upper = upper_wilson, precision = precision)
  }

  return(results)
}

# Apply BLB to each Region - FOR MEDIAN
all_blb_results <- list()
for (i in 1:nrow(empirical_data)) {
  all_blb_results[[empirical_data$Region[i]]] <- perform_blb_empirical(empirical_data[i, ], n_vals, B, b)
}

# Calculate Wilson CI and precision for each Region - FROM EMPIRICAL DATA
all_wilson_results <- list()
for (i in 1:nrow(empirical_data)) {
  all_wilson_results[[empirical_data$Region[i]]] <- calculate_wilson_ci_precision(empirical_data[i, ], n_vals)
}

# --- Data Preparation for Plotting ---
plot_data <- data.frame(
  n = numeric(),
  Prevalence = numeric(),
  Median = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  Precision = numeric(),
  Region = character(),
  stringsAsFactors = FALSE
)

for (region in empirical_data$Region) {
  for (n in n_vals) {
    n_char <- as.character(n)

    # Extract values
    median_val <- median(all_blb_results[[region]][[n_char]], na.rm = TRUE)
    ci_lower_val <- all_wilson_results[[region]][[n_char]]$ci_lower
    ci_upper_val <- all_wilson_results[[region]][[n_char]]$ci_upper
    precision_val <- all_wilson_results[[region]][[n_char]]$precision

    temp_data <- data.frame(
      n = n,
      Prevalence = empirical_data$proportion[empirical_data$Region == region] * 100,
      Median = median_val,
      CI_Lower = ci_lower_val,
      CI_Upper = ci_upper_val,
      Precision = precision_val,
      Region = region
    )

    plot_data <- bind_rows(plot_data, temp_data)
  }
}

# Remove rows with NA median (where n > dissected)
plot_data <- plot_data %>% filter(!is.na(Median))

# Reorder Region factor based on Prevalence
empirical_data <- empirical_data %>%
  arrange(desc(proportion))

# Convert Region to a factor for correct color mapping
plot_data$Region <- factor(plot_data$Region, levels = empirical_data$Region)

# Define a color palette
color_palette <- brewer.pal(n = 8, name = "Set1") # Increased palette size
names(color_palette) <- levels(plot_data$Region)

# Plot 1: Median prevalence and confidence intervals
plot_1 <- ggplot(plot_data, aes(x = n, y = Median, color = Region, group = Region)) +
      geom_hline(aes(yintercept = Prevalence, color = Region),
             linetype = "solid", linewidth = 0.5) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region, color = Region),
              alpha = 0.25,
              linetype = "dotted",
              linewidth = 0.2) +
   # Add vertical lines for n=25, 100 and 150
  annotate("segment", x = 20, xend = 20, y = 22, yend = 80,
           linewidth = 0.7, color = "black")  +
   annotate("segment", x = 30, xend = 30, y = 5, yend = 22,
           linewidth = 0.7, color = "black")+  
             geom_vline(xintercept = 110, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 160, color = "blue", linetype = "dashed") +
  annotate("rect", xmin = 111, xmax = 159, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.25)+
   scale_x_continuous(breaks = c(3,  20, 35, 50, 70, 90, 110, 125, 150, 180, 210, 250, 300, 350, 400),
                     labels = as.character(c(3,  20, 35, 50, 70, 90, 110, 125, 150, 180, 210, 250, 300, 350, 400)),
                     limits = c(0, 400)) +
 scale_y_continuous(breaks = c(0, 15, 25, 40, 50, 65, 80),
                     labels = as.character(c(0, 15, 25, 40, 50, 65, 80)))  +
 geom_point(size = 1) +
 geom_line(linetype = "dashed") +
 labs(
    x = "Sample size",
    y = "Prevalence (%)",
    tag = "D" # Add tag
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette)

# Plot 2: Precision
plot_2 <- ggplot(plot_data, aes(x = n, y = Precision, color = Region, group = Region)) +
  geom_line() +
  geom_point(size=0.9) +
  scale_y_continuous(breaks = c(4, 6, 10, 15, 20, 25, 30, 35, 40),
                     labels = as.character(c(4, 6, 10, 15, 20, 25, 30, 35, 40)),
                     limits = c(4, 40))  +
 scale_x_continuous(breaks = c(3,  20, 35, 50, 70, 90, 110, 125, 150, 180, 210, 250, 300, 350, 400),
                     labels = as.character(c(3,  20, 35, 50, 70, 90, 110, 125, 150, 180, 210, 250, 300, 350, 400)),
                     limits = c(0, 400)) +
                     labs(
    x = "Sample size",
    y = "Precision (%)",
    color = "Region",
    tag = "E"
  ) +
    theme_bw() +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1),
    legend.position = "none" # Hide legend on this plot
  ) +
  scale_color_manual(values = color_palette)

# ______Plot changes in precision variability as an outcome of the sample size increase

# Define the sample size ranges
ranges <- list(
 "3-9" = 3:9,
  "9-14" = 9:14,
  "14-19" = 14:19,
  "19-24" = 19:24,
  "24-30" = 24:30,
  "30-35" = 30:35,
  "35-40" = 35:40,
  "40-45" = 40:45,
"45-50" = 45:50,
"50-55" = 50:55,
"55-60" = 55:60,
"60-65" = 60:65,
"65-70" = 65:70,
"70-75" = 70:75,
"75-80" = 75:80,
"80-85" = 80:85,
"85-90" = 85:90,
"90-95" = 90:95,
"95-100" = 95:100,
"100-105" = 100:105,
"105-110" = 105:110,
"110-115" = 110:115,
"115-120" = 115:120,
"120-125" = 120:125,
"125-130" = 125:130,
"130-135" = 130:135,
"135-140" = 135:140,
"140-145" = 140:145,
"145-150" = 145:150,
"150-155" = 150:155,
"155-160" = 155:160,
"160-165" = 160:165,
"165-170" = 165:170,
"170-175" = 170:175,
"175-180" = 175:180,
"180-185" = 180:185,
"185-190" = 185:190,
"190-195" = 190:195,
"195-200" = 195:200,
"200-205" = 200:205,
"205-210" = 205:210,
"210-215" = 210:215,
"215-220" = 215:220,
"220-225" = 220:225,
"225-230" = 225:230,
"230-235" = 230:235,
"235-240" = 235:240,
"240-245" = 240:245,
"245-250" = 245:250,
"250-255" = 250:255,
"255-260" = 255:260,
"260-265" = 260:265,
"265-270" = 265:270
   )


# 1. Create a data frame with all combinations of Region and Range
all_combinations <- expand.grid(Region = unique(plot_data$Region), Range = names(ranges), stringsAsFactors = FALSE)

# 2. Create a lookup table to convert Range names to a vector of n values
range_lookup <- data.frame(Range = names(ranges), range_values = I(ranges))

# 3. Join the lookup table to all_combinations
all_combinations <- left_join(all_combinations, range_lookup, by = "Range")

# Function to calculate precision change, handling potential empty data
calculate_precision_change_safe <- function(region, range_values, plot_data) {
  region_data <- plot_data %>% filter(Region == region, n %in% range_values)
  if (nrow(region_data) < 2) {
    return(NA) # Or 0 if you want to assume zero change
  }
  return(max(region_data$Precision) - min(region_data$Precision))
}

# Calculate precision changes, using the safe version
precision_changes <- all_combinations %>%
  rowwise() %>%
  mutate(PrecisionChange = calculate_precision_change_safe(Region, range_values, plot_data)) %>%
  ungroup()

# Convert 'Range' to an ordered factor
precision_changes$Range <- factor(precision_changes$Range, levels = names(ranges))

# Create the plot 
plot_3<-ggplot(precision_changes, aes(x = Range, y = PrecisionChange, color = Region, group = Region)) +
  geom_line(na.rm = TRUE) +
  geom_point(size=1, na.rm = TRUE) +
  labs(
    x = "Sample size increment",
    y = "Change in precision",
    color = "Region",
    tag = "F"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  scale_color_manual(values = color_palette) +
  scale_x_discrete(breaks = levels(precision_changes$Range)[seq(1, length(levels(precision_changes$Range)), by = 5)])

# Combine plots horizontally
combined_plot <- plot_1 + plot_2 + plot_3 + plot_layout(guides = "collect")

# Output the combined plot
print(combined_plot)


#Export data in CSV
write.csv(plot_data, file = "plot_data.csv", row.names = T)

#Export data in CSV
write.csv(precision_changes, file = "precision_changes.csv", row.names = T)



