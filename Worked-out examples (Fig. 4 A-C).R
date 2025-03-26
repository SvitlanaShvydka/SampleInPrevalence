# Load necessary packages
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)



# Empirical prevalence data of acanthocephalan parasite Pomphorhynchus tereticollis (Rudolphi) from the definitive (cyprinid fish Phoxinus phoxinus (L.)) and the intermediate (crustacean amphipod Gammarus balcanicus Schäferna) host from the Malé Morské Oko lake, Eastern Slovakia, obtained in 2001 and 2022-2023 (S1 Supplementary data)
empirical_data <- data.frame(
  Host_and_Year = c("PP in 2001", "GB in 2001", "PP in 2022-2023", "GB in 2022-2023"),
  dissected = c(129, 612, 14, 630),
  infected = c(112, 280, 8, 25),
  proportion = c(0.87, 0.46, 0.57, 0.04)
)

# Parameters
set.seed(500910000)
B <- 1000  # Number of bootstrap samples
b <- 100   # Bootstrap resamples within each bootstrap sample
n_vals <- c(3, 6, 9, 12, 14, 19, 24, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300) # Subsample sizes

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

# Apply BLB to each host and year combination - FOR MEDIAN
all_blb_results <- list()
for (i in 1:nrow(empirical_data)) {
  all_blb_results[[empirical_data$Host_and_Year[i]]] <- perform_blb_empirical(empirical_data[i, ], n_vals, B, b)
}

# Calculate Wilson CI and precision for each host and year combination - FROM EMPIRICAL DATA
all_wilson_results <- list()
for (i in 1:nrow(empirical_data)) {
  all_wilson_results[[empirical_data$Host_and_Year[i]]] <- calculate_wilson_ci_precision(empirical_data[i, ], n_vals)
}

# --- Data Preparation for Plotting ---
plot_data <- data.frame(
  n = numeric(),
  Prevalence = numeric(),
  Median = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  Precision = numeric(),
  Host_and_Year = character(),
  stringsAsFactors = FALSE
)

for (host_year in empirical_data$Host_and_Year) {
  for (n in n_vals) {
    n_char <- as.character(n)

    # Extract values
    median_val <- median(all_blb_results[[host_year]][[n_char]], na.rm = TRUE)
    ci_lower_val <- all_wilson_results[[host_year]][[n_char]]$ci_lower
    ci_upper_val <- all_wilson_results[[host_year]][[n_char]]$ci_upper
    precision_val <- all_wilson_results[[host_year]][[n_char]]$precision

    temp_data <- data.frame(
      n = n,
      Prevalence = empirical_data$proportion[empirical_data$Host_and_Year == host_year] * 100,
      Median = median_val,
      CI_Lower = ci_lower_val,
      CI_Upper = ci_upper_val,
      Precision = precision_val,
      Host_and_Year = host_year
    )

    plot_data <- bind_rows(plot_data, temp_data)
  }
}

# Remove rows with NA median (where n > dissected)
plot_data <- plot_data %>% filter(!is.na(Median))

# Reorder Host_and_Year factor based on Prevalence
empirical_data <- empirical_data %>%
  arrange(desc(proportion))

# Convert Host_and_Year to a factor for correct color mapping
plot_data$Host_and_Year <- factor(plot_data$Host_and_Year, levels = empirical_data$Host_and_Year)

# Define a color palette
color_palette <- brewer.pal(n = 4, name = "Set1")
names(color_palette) <- levels(plot_data$Host_and_Year)

# Plot 1: Median prevalence and confidence intervals
plot_1 <- ggplot(plot_data, aes(x = n, y = Median, color = Host_and_Year, group = Host_and_Year)) +
  geom_hline(aes(yintercept = Prevalence, color = Host_and_Year),
             linetype = "solid", linewidth = 0.5) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper, fill = Host_and_Year, color = Host_and_Year),
              alpha = 0.25,
              linetype = "dotted",
              linewidth = 0.2) +
  # Add vertical lines for n=25, 100 and 150
  annotate("segment", x = 20, xend = 20, y = 25, yend = 67,
           linewidth = 0.7, color = "black")  +
   annotate("segment", x = 100, xend = 100, y = 0.5, yend = 9.5,
           linewidth = 0.7, color = "black")+  
  annotate("segment", x = 30, xend = 30, y = 70, yend = 95,
           linewidth = 0.7, color = "black")+
           geom_vline(xintercept = 110, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 135, color = "blue", linetype = "dashed") +
  annotate("rect", xmin = 111, xmax = 134, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.25) +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  scale_x_continuous(breaks = c(3, 20, 40, 60, 80, 100, 125, 150, 180),
                     labels = as.character(c(3, 20, 40, 60, 80, 100, 125, 150, 180)),
                     limits = c(0, 180)) +
  scale_y_continuous(breaks = c(5, 15, 25, 35,  45,  55,  65,  75, 85,  95),
                     labels = as.character(c(5, 15, 25, 35,  45,  55,  65,  75, 85,  95)),
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

# -Plot precision
plot_2 <- ggplot(plot_data, aes(x = n, y = Precision, color = Host_and_Year, group = Host_and_Year)) +
  geom_line() +  # Use linewidth instead of size for smooth curves
  geom_point(size=0.9) +  # Keep using size for points
  scale_color_manual(values = color_palette) +  # Ensure consistent colors
  scale_x_continuous(breaks = c(3, 20, 40, 60, 80, 100, 125, 150, 200, 250, 300),
                     labels = as.character(c(3, 20, 40, 60, 80, 100, 125, 150, 200, 250, 300)),
                     limits = c(3, 300))  +
  scale_y_continuous(breaks = c(2, 4, 6, 10, 15, 20, 25, 30, 35, 40),
                     labels = as.character(c(2, 4, 6, 10, 15, 20, 25, 30, 35, 40)),
                     limits = c(2, 40))  +
  labs(
    x = "Sample size",
    y = "Precision (%)",
    color = "Host and Year",
    tag = "B"  # Add tag for plot identification
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels if needed
    legend.position = "false"  # Adjust legend position
  )


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
# (3, 6, 9, 12, 14, 17, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 115, 130, 150, 180, 210, 240, 270, 300)
# Function to calculate precision change within a range
calculate_precision_change <- function(data, range_values) { # Simplified function
  if (nrow(data) < 2) {
    return(NA)  # Not enough data points to calculate change
  }
  precision_change <- max(data$Precision) - min(data$Precision)
  return(precision_change)
}

# Calculate precision changes for each range and Host_and_Year
precision_changes <- do.call(rbind, lapply(unique(plot_data$Host_and_Year), function(region) { # Apply lapply and rbind
  region_data <- plot_data %>% filter(Host_and_Year == region)
  data.frame(
    Host_and_Year = region,
    Range = names(ranges),
    PrecisionChange = sapply(ranges, function(range_values) calculate_precision_change(filter(region_data, n %in% range_values)))
  )
}))


# Convert 'Range' to an ordered factor with strict ordering
precision_changes$Range <- factor(precision_changes$Range, levels = names(ranges))

# Create the line plot (curves)
plot_3 <- ggplot(precision_changes, aes(x = Range, y = PrecisionChange, color = Host_and_Year, group = Host_and_Year)) +
  geom_line() +  # Add lines to connect the points
  geom_point(size=1) + # Add points at each data point
 scale_x_discrete(breaks = levels(precision_changes$Range)[seq(1, length(levels(precision_changes$Range)), by = 5)])+
 labs(
    x = "Sample size increment",
    y = "Change in precision",
    color = "Host_and_Year",
    tag = "C" # Add tag
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = color_palette)



# Combine graphs horizontally
combined_plot <- plot_1 + plot_2 + plot_3 + plot_layout(guides = "collect")

# Output the combined plot
print(combined_plot)


#Export data in CSV
write.csv(plot_data, file = "plot_data_precision-Teraticol.csv", row.names = T)


#Export data in CSV
write.csv(precision_changes, file = "precision_changes_terat.csv", row.names = T)

