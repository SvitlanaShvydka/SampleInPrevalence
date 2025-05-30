# Load necessary packages
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(tidyr)

# Empirical data on the prevalence of Echinococcus multilocularis in red foxes from eight regions in Slovakia were accumulated over 5 years (2000, 2001, 2002, 2005/2006, and 2007), as reported by Miterpáková and Dubinský (2011) DOI: 10.2478/s11687-011-0023-5
empirical_data <- data.frame(
  Region = c("ZA", "TN", "PO", "BB", "NR", "KE", "TT", "BA"),
  dissected = c(379, 254, 647, 762, 476, 670, 271, 172),
  infected = c(197, 105, 259, 228, 118, 127, 42, 22),
  proportion = c(0.520, 0.413, 0.400, 0.299, 0.248, 0.190, 0.155, 0.128),
  Population = c(6809, 4502, 8973, 9454, 6344, 6754, 4146, 2053)
)

# Parameters
n_vals <- c(3, 6, 9, 12, 14, 19, 24, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300, 305, 310, 315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390, 395, 400, 405, 410, 415, 420)

# Population sizes
population_sizes <- c("BA" = 1481, "TT" = 2991, "NR" = 4577, "TN" = 3248, "ZA" = 4913, "BB" = 6821, "KE" = 4873, "PO" = 6474)

# CONF.prop function (Wilson CI with population correction)
CONF.prop <- function(alpha, x = NULL,
                      sample.prop = mean(x), n = length(x),
                      N = Inf, unsampled = FALSE) {
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1) {
    stop('Error: alpha should be a single numeric value between 0 and 1')
  }
  if ((!missing(x) && !missing(sample.prop)) || (missing(x) && missing(sample.prop))) {
    stop('Error: specify data or sample proportion but not both')
  }
  if (!missing(x)) {
    if (is.logical(x)) x <- x + 0
    if (!is.numeric(x) || !all(x %in% c(0L, 1L))) {
      stop('Error: x should be binary data')
    }
    P <- mean(x)
    n <- length(x)
  } else {
    P <- sample.prop
  }
  if (!is.numeric(P) || length(P) != 1) {
    stop('Error: sample proportion should be a single numeric value')
  }
  if (!is.numeric(n) || as.integer(n) != n || length(n) != 1 || n < 1) {
    return(list(lower = NA, upper = NA))
  }
  if (alpha == 0) return(list(lower = 0, upper = 1))
  if (alpha == 1) return(list(lower = P, upper = P))

  if ((alpha > 0) && (alpha < 1)) {
    if (N == Inf) {
      psi <- qchisq(1 - alpha, df = 1) / (2 * n)
    } else {
      if (unsampled) {
        psi <- ((N - 1) / (N - n)) * qchisq(1 - alpha, df = 1) / (2 * n)
      } else {
        psi <- ((N - n) / (N - 1)) * qchisq(1 - alpha, df = 1) / (2 * n)
      }
    }
    T1 <- (P + psi) / (1 + 2 * psi)
    T2 <- sqrt(2 * psi * P * (1 - P) + psi^2) / (1 + 2 * psi)
    L <- T1 - T2
    U <- T1 + T2
    return(list(lower = L, upper = U))
  } else {
    return(list(lower = NA, upper = NA))
  }
}

# Function to calculate Wilson CI and precision - FROM EMPIRICAL DATA
calculate_wilson_ci_precision <- function(data_row, n_vals, pop_size = Inf) {
  results <- list()

  for (n in n_vals) {
    if (n > data_row$dissected) {
      results[[as.character(n)]] <- list(ci_lower = NA, ci_upper = NA, precision = NA)
      next
    }

    # Calculate Wilson score interval using CONF.prop
    wilson_ci <- CONF.prop(
      alpha = 0.05,
      sample.prop = data_row$infected / data_row$dissected,
      n = n,
      N = pop_size  # Use the provided population size
    )

    lower_wilson <- wilson_ci$lower * 100
    upper_wilson <- wilson_ci$upper * 100

    # Calculate precision (width of the CI)
    precision <- (upper_wilson - lower_wilson) / 2

    results[[as.character(n)]] <- list(ci_lower = lower_wilson, ci_upper = upper_wilson, precision = precision)
  }

  return(results)
}

# Calculate Wilson CI and precision for each host and year combination
all_wilson_results <- list()
for (i in 1:nrow(empirical_data)) {
  # Determine the population size for the current Host_and_Year
  pop_size <- ifelse(empirical_data$Region[i] %in% names(population_sizes),
                     population_sizes[empirical_data$Region[i]],
                     Inf)

  all_wilson_results[[empirical_data$Region[i]]] <- calculate_wilson_ci_precision(empirical_data[i, ], n_vals, pop_size)
}

# --- Data Preparation for Plotting ---
plot_data <- data.frame(
  n = numeric(),
  Prevalence = numeric(),
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
    ci_lower_val <- all_wilson_results[[region]][[n_char]]$ci_lower
    ci_upper_val <- all_wilson_results[[region]][[n_char]]$ci_upper
    precision_val <- all_wilson_results[[region]][[n_char]]$precision

    temp_data <- data.frame(
      n = n,
      Prevalence = empirical_data$proportion[empirical_data$Region == region] * 100,
      CI_Lower = ci_lower_val,
      CI_Upper = ci_upper_val,
      Precision = precision_val,
      Region = region
    )

    plot_data <- bind_rows(plot_data, temp_data)
  }
}

# Print the final plot_data
print(plot_data)

#Export data in CSV
write.csv(plot_data, file = "Adjusted CI-Echinoc.csv", row.names = T)

# Reorder Host_and_Year factor based on Prevalence
empirical_data <- empirical_data %>%
  arrange(desc(proportion))

# Convert Host_and_Year to a factor for correct color mapping
plot_data$Region <- factor(plot_data$Region, levels = empirical_data$Region)

# Define a color palette
color_palette <- RColorBrewer::brewer.pal(n = 8, name = "Set1")
names(color_palette) <- levels(plot_data$Region)

# Plot 1: Prevalence and confidence intervals
ggplot(plot_data, aes(x = n, color = Region, group = Region)) +
  geom_hline(aes(yintercept = Prevalence, color = Region),
             linetype = "solid", linewidth = 0.5) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region, color = Region),
              alpha = 0.25,
              linetype = "dotted",
              linewidth = 0.2) +
  # Add vertical lines for n=20, 30 and 100
  annotate("segment", x = 145, xend = 175, y = 32.4, yend = 60,
           linewidth = 0.7, color = "black", arrow = arrow(length = unit(0.3, "cm"), type = "open", ends = "first"))  +
     annotate("segment", x = 330, xend = 360, y = 35, yend = 60,
           linewidth = 0.7, color = "black", arrow = arrow(length = unit(0.3, "cm"), type = "open", ends = "first"))+
    scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  scale_x_continuous(breaks = c(5, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420),
                     labels = as.character(c(5, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420)),
                     limits = c(0, 420)) +
  scale_y_continuous(breaks = c(5, 15, 25, 35,  45,  55,  65,  75, 85),
                     labels = as.character(c(5, 15, 25, 35,  45,  55,  65,  75, 85)),
                     limits = c(0, 85))+
      labs(
    x = "Sample size",
    y = "Prevalence (%)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
