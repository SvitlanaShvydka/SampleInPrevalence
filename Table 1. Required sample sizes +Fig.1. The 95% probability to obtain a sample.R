# Load necessary packages
library(binom)  
library(dplyr)  
library(tibble) 
library(ggplot2)
library(RColorBrewer)
library(tidyr)

# Set prevalence values (as percentages)
prevalence_values <- c(99, 98, 97, 96, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 4, 3, 2, 1)

# Create an empty data frame to store the results
results_df <- data.frame(
  Prevalence = numeric(),
  SampleSize_D_leq_P_or_100_minus_P = numeric(),
  Precision_D_leq_P_or_100_minus_P = numeric(),
  SampleSize_nP_and_n100_P_geq5 = numeric(),
  Precision_nP_and_n100_P_geq5 = numeric(),
  Precision_nP_and_n100_P_geq5_Wald = numeric(),
  SampleSize_P_neq0_and_P_neq100 = numeric(),
  Precision_P_neq0_and_P_neq100 = numeric(),
  Target_precision = numeric(),
  stringsAsFactors = FALSE
)

# Сalculate Wilson score interval
wilson_interval_calc <- function(n, p_percent, confidence_level = 0.95) {
  # Convert percentage to proportion for calculation
  p <- p_percent / 100
  if (is.na(n) || n <= 0) {
    return(list(lower = NA, upper = NA))
  }
  x <- round(n * p, 2)
  wilson_interval <- tryCatch({
    binom.confint(x, n, conf.level = confidence_level, methods = "wilson")
  }, error = function(e) {
    return(data.frame(lower = NA, upper = NA))
  })
  if (all(is.na(wilson_interval))) {
    return(list(lower = NA, upper = NA))
  }
  return(list(lower = wilson_interval[1, "lower"], upper = wilson_interval[1, "upper"]))
}

# Calculate margin of error for Wilson score interval 
wilson_precision <- function(lower, upper) {
  if (is.na(lower) || is.na(upper)) {
    return(NA)
  }
  return((upper - lower) / 2)
}

# Calculate margin of error for Wald interval 
wald_precision <- function(n, p_percent) {
  p <- p_percent / 100
  if (is.na(n) || n <= 0) {
    return(NA)
  }
  return(1.96 * sqrt((p * (1 - p)) / n) * 100)
}

# Loop through the prevalence values
for (p_percent in prevalence_values) {
  # Convert percentage to proportion for calculation
  p <- p_percent / 100

  # Calculate target precision based on the condition
  if (p_percent < 50) {
    target_precision <- p_percent
    target_prevalence <- p
  } else {
    target_precision <- 100 - p_percent
    target_prevalence <- target_precision / 100
  }

  # -----------------------------------------------------------------------
  # --- 1. Calculate n where D≤P or D≤100-P with 95% probability (Wald + Wilson Adj CIs) ---
  # -----------------------------------------------------------------------

  # --- 1a. Estimate n using Wald margin of error ---
  wald_n <- tryCatch({
    if (p_percent < 50) {
      tibble::tibble(
        n = seq_len(1e6),
        margin_error = 1.96 * sqrt((p * (1 - p)) / n)
      ) |>
        dplyr::filter(margin_error <= p) |>
        dplyr::pull(n) %>%
        .[1]
    } else {
      tibble::tibble(
        n = seq_len(1e6),
        margin_error = 1.96 * sqrt((p * (1 - p)) / n)
      ) |>
        dplyr::filter(margin_error <= (1 - p)) |>
        dplyr::pull(n) %>%
        .[1]
    }
  }, error = function(e) {
    NA
  })

  # --- 1b. Iteratively increase n until Wilson precision is met ---
  n <- wald_n  # Start with the Wald estimate
  precision_met <- FALSE
  max_iterations <- 1000
  iteration <- 0

  while (!precision_met && iteration < max_iterations) {
    iteration <- iteration + 1
    wilson_interval_result <- wilson_interval_calc(n, p_percent)
    precision <- ifelse(is.na(n), NA, wilson_precision(wilson_interval_result$lower, wilson_interval_result$upper) * 100)

    if (!is.na(precision) && precision <= target_precision) {
      precision_met <- TRUE
    } else {
      n <- n + 10
    }

    if (n > 1e6) {
      n <- NA
      precision <- NA
      break
    }
  }

  precision_d_leq_p <- precision # The precision after Wilson adjustment

  # -----------------------------------------------------------------------
  # --- 2. Calculate n where nP >= 5 and n(1-P) >= 5 with 95% probability ---
  # -----------------------------------------------------------------------

  if (p_percent < 50) {
    result_np5 <- tryCatch({
      tibble::tibble(n = seq_len(1e6), pr = stats::pbinom(4, n, p, lower.tail = FALSE)) |>
        dplyr::filter(pr >= 0.95) |>
        dplyr::pull(n) %>%
        .[1]
    }, error = function(e) {
      NA
    })
  } else {
    result_np5 <- tryCatch({
      tibble::tibble(n = seq_len(1e6), pr = stats::pbinom(4, n, target_prevalence, lower.tail = FALSE)) |>
        dplyr::filter(pr >= 0.95) |>
        dplyr::pull(n) %>%
        .[1]
    }, error = function(e) {
      NA
    })
  }

  # Calculate precision for result_np5 using Wilson and Wald CIs
  wilson_interval_np5 <- wilson_interval_calc(result_np5, p_percent)
  precision_np5 <- ifelse(is.na(result_np5), NA, wilson_precision(wilson_interval_np5$lower, wilson_interval_np5$upper) * 100)
  precision_np5_wald <- ifelse(is.na(result_np5), NA, wald_precision(result_np5, p_percent))

  # -----------------------------------------------------------------------
  # --- 3. Calculate n where P != 0% and P != 100% with 95% probability ---
  # -----------------------------------------------------------------------

  result_pneq0 <- NA
  if (p > 0 && p < 1) {
    if (p_percent < 50) {
      result_pneq0 <- tryCatch({
        tibble::tibble(n = seq_len(1e6), pr = stats::pbinom(0, n, p, lower.tail = FALSE)) |>
          dplyr::filter(pr >= 0.95) |>
          dplyr::pull(n) %>%
          .[1]
      }, error = function(e) {
        NA
      })
    } else {
      result_pneq0 <- tryCatch({
        tibble::tibble(n = seq_len(1e6), pr = stats::pbinom(0, n, target_prevalence, lower.tail = FALSE)) |>
          dplyr::filter(pr >= 0.95) |>
          dplyr::pull(n) %>%
          .[1]
      }, error = function(e) {
        NA
      })
    }
  }

  # Calculate precision for result_pneq0 using Wilson CIs
  wilson_interval_pneq0 <- wilson_interval_calc(result_pneq0, p_percent)
  precision_pneq0 <- ifelse(is.na(result_pneq0), NA, wilson_precision(wilson_interval_pneq0$lower, wilson_interval_pneq0$upper) * 100)

  # Store the results in the data frame
  new_row <- data.frame(
    Prevalence = p_percent,
    SampleSize_D_leq_P_or_100_minus_P = n,
    Precision_D_leq_P_or_100_minus_P = round(precision_d_leq_p, 1),
    SampleSize_nP_and_n100_P_geq5 = result_np5,
    Precision_nP_and_n100_P_geq5 = round(precision_np5, 1),
    Precision_nP_and_n100_P_geq5_Wald = round(precision_np5_wald, 1),
    SampleSize_P_neq0_and_P_neq100 = result_pneq0,
    Precision_P_neq0_and_P_neq100 = round(precision_pneq0, 1),
    Target_precision = target_precision
  )

  results_df <- rbind(results_df, new_row)
}


#Export data in CSV
write.csv(results_df, file = "n+precis.csv", row.names = T)

#-------------------finite populations adjusted n
# Define vector of finite population sizes
population_sizes <- c(100, 500, 1000, 5000, 10000)

# Create a data frame to store the adjusted sample sizes
adjusted_sample_sizes_df <- data.frame(
  Population_Size = numeric(),
  Prevalence = numeric(),
  Adjusted_SampleSize_D_leq_P_or_100_minus_P = numeric(),
  Adjusted_SampleSize_nP_and_n100_P_geq5 = numeric(),
  Adjusted_SampleSize_P_neq0_and_P_neq100 = numeric(),
  Target_precision = numeric(),
  stringsAsFactors = FALSE
)

# Create a data frame to store the adjusted Wilson CIs and adjusted n
adjusted_results_df <- data.frame(
  Population_Size = numeric(),
  Prevalence = numeric(),
  Adjusted_SampleSize_D_leq_P_or_100_minus_P = numeric(),
  Wilson_Precision_D_leq_P_or_100_minus_P = numeric(),
  Adjusted_SampleSize_nP_and_n100_P_geq5 = numeric(),
  Wilson_Precision_nP_and_n100_P_geq5 = numeric(),
  Adjusted_SampleSize_P_neq0_and_P_neq100 = numeric(),
  Wilson_Precision_P_neq0_and_P_neq100 = numeric(),
  Target_precision = numeric(),
  stringsAsFactors = FALSE
)


# CONF.prop function (Wilson score interval)
#The idea from O’NEILL, B. AND FULTZ, N. (2020) stat.extend: highest density regions and other functions of distributions. R package, Version 0.1.4. https://CRAN.R-project.org/package=stat.extend

CONF.prop <- function(alpha, x = NULL, 
                      sample.prop = mean(x), n = length(x),
                      N = Inf, unsampled = FALSE) {
  
  #Check input alpha
  if (!is.numeric(alpha))   { stop('Error: alpha should be numeric') }
  if (length(alpha) != 1)   { stop('Error: alpha should be a single value'); }
  if (alpha < 0)            { stop('Error: alpha is negative'); }
  if (alpha > 1)            { stop('Error: alpha is greater than one'); }
  
  #Check congruence of data inputs
  if ((!missing(x) && !missing(sample.prop))) {
    if (abs(sample.prop - mean(x)) < 1e-15) 
      warning('specify data or sample proportion but not both') else
        stop('Error: specify data or sample proportion but not both'); }
  if ((!missing(x) && !missing(n))) {
    if (n != length(x)) 
      stop('Error: specify data or sample proportion but not both'); }
  
  #Check data inputs
  P <- sample.prop;
  if (!missing(x)) {
    xexpr <- deparse(substitute(x));
    if (is.logical(x)) x <- x + 0;
    if (!is.numeric(x))     { stop('Error: x should be numeric') }
    if (!all(x %in% c(0L,1L))) { stop('Error: x should be binary data') } }
  if (!is.numeric(P))       { stop('Error: sample proportion should be numeric') }
  if (length(P) != 1)       { stop('Error: sample proportion should be a single value'); }
  if (!is.numeric(n))       { stop('Error: n should be numeric') }
  if (as.integer(n) != n)   { stop('Error: n should be an integer') }
  if (length(n) != 1)       { stop('Error: n should be a single value'); }
  if (n < 1)                { return(list(lower = NA, upper = NA)) } # Modified this line
  
  #Compute the confidence interval in trivial cases
  if (alpha == 0) {
    return(list(lower = 0, upper = 1)) } # Modified this line
  if (alpha == 1) {
    return(list(lower = P, upper = P)) } # Modified this line
  
  #Compute the confidence interval in non-trivial case
  if ((alpha > 0) && (alpha < 1)) {
    if (N == Inf) {
      psi <- qchisq(1-alpha, df = 1)/(2*n) } else {
        if (unsampled) {
          psi <- ((N-1)/(N-n))*qchisq(1-alpha, df = 1)/(2*n) } else {
            psi <- ((N-n)/(N-1))*qchisq(1-alpha, df = 1)/(2*n) } }
    T1  <- (P+psi)/(1+2*psi);
    T2  <- sqrt(2*psi*P*(1-P) + psi^2)/(1+2*psi);
    L   <- T1 - T2;
    U   <- T1 + T2;
    return(list(lower = L, upper = U))
  } else {
    return(list(lower = NA, upper = NA)) # Added this line
  }
}

# Function to adjust the sample size based on the population size, accounting for finite population correction
adjust_sample_size <- function(n, population_size) {
  if (!is.na(n)) { # Modified this line
    n <- as.numeric(n) # Added this line
    if (n > 0 && n <= population_size) {
      # Use the specified formula for adjusted sample size
      adjusted_n <- ceiling(n / (1 + n / population_size))
      return(adjusted_n)
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}

# Function to calculate precision from Wilson CI
calculate_precision <- function(lower, upper) {
  if (is.na(lower) || is.na(upper)) {
    return(NA)
  }
  return(((upper - lower) / 2) * 100) # Return precision as percentage
}

# Loop through the population sizes
for (pop_size in population_sizes) {
  # Loop through the rows of the original results data frame
  for (i in 1:nrow(results_df)) {
    p_percent <- results_df$Prevalence[i] / 100  # Convert percentage to proportion
    n_d_leq_p <- results_df$SampleSize_D_leq_P_or_100_minus_P[i]
    n_np5 <- results_df$SampleSize_nP_and_n100_P_geq5[i]
    n_pneq0 <- results_df$SampleSize_P_neq0_and_P_neq100[i]
    target_precision <- results_df$Target_precision[i]

    # Adjust n for finite population
    adjusted_n_d_leq_p <- adjust_sample_size(n_d_leq_p, pop_size)
    adjusted_n_np5 <- adjust_sample_size(n_np5, pop_size)
    adjusted_n_pneq0 <- adjust_sample_size(n_pneq0, pop_size)

    # Store adjusted sample sizes in the first data frame
    new_sample_sizes_row <- data.frame(
      Population_Size = pop_size,
      Prevalence = p_percent * 100,
      Adjusted_SampleSize_D_leq_P_or_100_minus_P = adjusted_n_d_leq_p,
      Adjusted_SampleSize_nP_and_n100_P_geq5 = adjusted_n_np5,
      Adjusted_SampleSize_P_neq0_and_P_neq100 = adjusted_n_pneq0,
      Target_precision = target_precision
    )
    adjusted_sample_sizes_df <- rbind(adjusted_sample_sizes_df, new_sample_sizes_row)

    # Calculate Wilson CIs
    wilson_d_leq_p <- tryCatch({
      CONF.prop(alpha = 0.05, sample.prop = p_percent, n = adjusted_n_d_leq_p, N = pop_size)
    }, error = function(e) {
      list(lower = NA, upper = NA)
    })
    wilson_np5 <- tryCatch({
      CONF.prop(alpha = 0.05, sample.prop = p_percent, n = adjusted_n_np5, N = pop_size)
    }, error = function(e) {
      list(lower = NA, upper = NA)
    })
    wilson_pneq0 <- tryCatch({
      CONF.prop(alpha = 0.05, sample.prop = p_percent, n = adjusted_n_pneq0, N = pop_size)
    }, error = function(e) {
      list(lower = NA, upper = NA)
    })

    # Calculate precision
    precision_d_leq_p <- calculate_precision(wilson_d_leq_p$lower, wilson_d_leq_p$upper)
    precision_np5 <- calculate_precision(wilson_np5$lower, wilson_np5$upper)
    precision_pneq0 <- calculate_precision(wilson_pneq0$lower, wilson_pneq0$upper)

    # Store adjusted n and precision in the second data frame
    new_results_row <- data.frame(
      Population_Size = pop_size,
      Prevalence = p_percent * 100,
      Adjusted_SampleSize_D_leq_P_or_100_minus_P = adjusted_n_d_leq_p,
      Wilson_Precision_D_leq_P_or_100_minus_P = precision_d_leq_p,
      Adjusted_SampleSize_nP_and_n100_P_geq5 = adjusted_n_np5,
      Wilson_Precision_nP_and_n100_P_geq5 = precision_np5,
      Adjusted_SampleSize_P_neq0_and_P_neq100 = adjusted_n_pneq0,
      Wilson_Precision_P_neq0_and_P_neq100 = precision_pneq0,
      Target_precision = target_precision
    )
    adjusted_results_df <- rbind(adjusted_results_df, new_results_row)
  }
}

#---------Prepare data for plotting


# Filter for Population_Size = 100
df_100 <- adjusted_results_df %>%
  filter(Population_Size == 100) %>%
  select(Prevalence, Adjusted_SampleSize_D_leq_P_or_100_minus_P) %>%
  rename(SampleSize_100 = Adjusted_SampleSize_D_leq_P_or_100_minus_P)

# Filter for Population_Size = 500
df_500 <- adjusted_results_df %>%
  filter(Population_Size == 500) %>%
  select(Prevalence, Adjusted_SampleSize_D_leq_P_or_100_minus_P) %>%
  rename(SampleSize_500 = Adjusted_SampleSize_D_leq_P_or_100_minus_P)

# Filter for Population_Size = 1000
df_1000 <- adjusted_results_df %>%
  filter(Population_Size == 1000) %>%
  select(Prevalence, Adjusted_SampleSize_D_leq_P_or_100_minus_P) %>%
  rename(SampleSize_1000 = Adjusted_SampleSize_D_leq_P_or_100_minus_P)

# Filter for Population_Size = 5000
df_5000 <- adjusted_results_df %>%
  filter(Population_Size == 5000) %>%
  select(Prevalence, Adjusted_SampleSize_D_leq_P_or_100_minus_P) %>%
  rename(SampleSize_5000 = Adjusted_SampleSize_D_leq_P_or_100_minus_P)

# Filter for Population_Size = 10000
df_10000 <- adjusted_results_df %>%
  filter(Population_Size == 10000) %>%
  select(Prevalence, Adjusted_SampleSize_D_leq_P_or_100_minus_P) %>%
  rename(SampleSize_10000 = Adjusted_SampleSize_D_leq_P_or_100_minus_P)

# Merge all data frames based on Prevalence
merged_df <- df_100 %>%
  full_join(df_500, by = "Prevalence") %>%
  full_join(df_1000, by = "Prevalence") %>%
  full_join(df_5000, by = "Prevalence") %>%
  full_join(df_10000, by = "Prevalence")



#---------------Plotting


# Define color palette
color_palette <- brewer.pal(n = 3, name = "Set1")

# Prepare original data for plotting
plotting_df <- results_df %>%
  dplyr::select(
    Prevalence,
    SampleSize_D_leq_P_or_100_minus_P,
    SampleSize_nP_and_n100_P_geq5,
    SampleSize_P_neq0_and_P_neq100
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("SampleSize"),
    names_to = "Method",
    values_to = "SampleSize"
  )

# Prepare merged_df for plotting
plotting_df_merged <- merged_df %>%
  tidyr::pivot_longer(
    cols = starts_with("SampleSize_"),
    names_to = "Method",
    values_to = "SampleSize"
  ) %>%
  mutate(Method = recode(Method,
                         "SampleSize_100" = "D≤P or D≤100-P, N=100",
                         "SampleSize_500" = "D≤P or D≤100-P, N=500",
                         "SampleSize_1000" = "D≤P or D≤100-P, N=1000",
                         "SampleSize_5000" = "D≤P or D≤100-P, N=5000",
                         "SampleSize_10000" = "D≤P or D≤100-P, N=10000")) %>%
  filter(!Method %in% c("D≤P or D≤100-P, N=500", "D≤P or D≤100-P, N=5000"))

# Combine original and merged data frames
plotting_df <- bind_rows(plotting_df, plotting_df_merged)

# Define labels for the legend
labels <- c(
  "SampleSize_nP_and_n100_P_geq5" = "x ≥ 5 and n - x ≥ 5",
  "SampleSize_P_neq0_and_P_neq100" = "P≠0%_and_P≠100%",
  "SampleSize_D_leq_P_or_100_minus_P" = "D≤P or D≤100-P",
  "D≤P or D≤100-P, N=100" = "D≤P or D≤100-P, N=100",
  "D≤P or D≤100-P, N=1000" = "D≤P or D≤100-P, N=1000",
  "D≤P or D≤100-P, N=10000" = "D≤P or D≤100-P, N=10000"
)

# Define line styles for the plot
line_styles <- c(
  "SampleSize_nP_and_n100_P_geq5" = "solid",
  "SampleSize_P_neq0_and_P_neq100" = "dotted",
  "SampleSize_D_leq_P_or_100_minus_P" = "dashed",
  "D≤P or D≤100-P, N=100" = "solid",
  "D≤P or D≤100-P, N=1000" = "solid",
  "D≤P or D≤100-P, N=10000" = "solid"
)

# Define colors for the plot
colors <- c(
  "SampleSize_nP_and_n100_P_geq5" = "black",
  "SampleSize_P_neq0_and_P_neq100" = "black",
  "SampleSize_D_leq_P_or_100_minus_P" = "black",
  "D≤P or D≤100-P, N=100" = color_palette[1],
  "D≤P or D≤100-P, N=1000" = color_palette[2],
  "D≤P or D≤100-P, N=10000" = color_palette[3]
)

# Create the plot
ggplot(plotting_df, aes(x = Prevalence, y = SampleSize, color = Method, linetype = Method, group = Method)) +
  geom_line(na.rm = TRUE) +
  geom_point(size = 0.8, na.rm = TRUE) +
  labs(
    x = "Prevalence (%)",
    y = "Sample size",
    color = "Method for estimating n",
    linetype = "Method for estimating n"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  scale_color_manual(values = colors, labels = labels) +
  scale_linetype_manual(values = line_styles, labels = labels) +
  scale_x_continuous(
    breaks = c(1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99),
    labels = as.character(c(1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99))
  ) +
scale_y_continuous(
    breaks = seq(0, 900, by = 30),
    limits = c(4, 920))   + 
    theme(legend.position = c(0.5, 0.72)) # Place the legend inside the plot, top center
  

#Export data in CSV
write.csv(adjusted_results_df, file = "n-adjasted.csv", row.names = T)

