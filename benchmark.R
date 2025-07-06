# Compare lambda.1se rule behavior for Cox regressions in glmnet 4.1-8 vs. 4.1-9
# Examine the effect of CV error normalization changes on model selection

library(doFuture)
library(ggplot2)

# Srcery palette
pal_srcery <- c("#EF2F27", "#519F50", "#FBB829", "#2C78BF", "#E02C6D", "#0AAEB3")

# Simulate one dataset and run regularized Cox regression ----
#' @param n Number of observations.
#' @param p Number of variables.
#' @param nzc Number of non-zero coefficients.
#' @param prob Censoring probability.
#' @param alpha Elastic-net mixing parameter.
#' @param use_1se Use 1 standard error rule?
#'
#' @return List of metrics.
simulate_cox <- function(
    n = 1000, p = 100, nzc = 10,
    prob, alpha,
    use_1se = TRUE) {
  # Example derived from ?glmnet::cv.glmnet
  x <- matrix(rnorm(n * p), n, p)

  beta <- rnorm(nzc)
  fx <- x[, seq(nzc)] %*% beta / 3
  hx <- exp(fx)
  ty <- rexp(n, hx)

  tcens <- rbinom(n = n, prob = prob, size = 1)
  y <- cbind(time = ty, status = 1 - tcens)

  fit_cv <- glmnet::cv.glmnet(x, y, family = "cox", nfolds = 5, alpha = alpha)
  lambda_selected <- if (use_1se) fit_cv$lambda.1se else fit_cv$lambda.min
  fit <- glmnet::glmnet(x, y, family = "cox", alpha = alpha, lambda = lambda_selected)

  selected_vars <- which(as.vector(fit$beta) != 0)
  true_vars <- seq(nzc)

  metrics <- list(
    total_selected = length(selected_vars),
    true_positives = sum(selected_vars %in% true_vars),
    false_positives = sum(!(selected_vars %in% true_vars)),
    prob = prob,
    alpha = alpha,
    lambda_type = if (use_1se) "lambda.1se" else "lambda.min"
  )

  metrics
}

# Install a specific glmnet version at runtime ----
install_glmnet <- function(version) {
  glmnet_installed <- installed.packages()[, "Package"] == "glmnet"
  if (any(glmnet_installed)) remove.packages("glmnet")
  message(paste("Installing glmnet", version), "...")
  remotes::install_version("glmnet", version = version, quiet = TRUE, upgrade = "never", force = TRUE)
}

# Run simulation on a specific glmnet version in parallel ----
#' @param version glmnet version to use.
#' @param n_reps Number of repetitions.
#'
#' @return Results data frame.
simulate_version <- function(version = c("4.1-8", "4.1-9"), n_reps) {
  install_glmnet(version)
  message(paste("Running experiments with glmnet version:", packageVersion("glmnet")))

  # Parameter grid
  prob_values <- seq(0.1, 0.9, by = 0.1)
  alpha_values <- c(1, 0.5)

  plan(multisession, workers = parallelly::availableCores() - 1)

  set.seed(42)
  results_list <- foreach(
    rep = 1:n_reps,
    .options.future = list(seed = TRUE)
  ) %dofuture% {
    results_single <- data.frame()
    for (prob in prob_values) {
      for (alpha in alpha_values) {
        metrics <- simulate_cox(prob = prob, alpha = alpha, use_1se = TRUE)
        results_single <- rbind(
          results_single,
          data.frame(
            rep = rep,
            prob = metrics$prob,
            alpha = metrics$alpha,
            total_selected = metrics$total_selected,
            true_positives = metrics$true_positives,
            false_positives = metrics$false_positives,
            version = version
          )
        )
      }
    }
    results_single
  }

  plan(sequential)

  results <- dplyr::bind_rows(results_list)

  as.data.frame(results)
}

# Ridgeline plots showing selection metric densities ----
plot_densities <- function(results_df) {
  # Convert to factor for proper ordering
  results_df$prob_factor <- factor(
    results_df$prob,
    levels = seq(0.1, 0.9, by = 0.1),
    ordered = TRUE
  )

  # Create plot for each version and alpha combination
  plots_list <- list()

  for (version in unique(results_df$version)) {
    for (alpha in unique(results_df$alpha)) {
      plot_data <- results_df |>
        dplyr::filter(version == !!version, alpha == !!alpha) |>
        tidyr::pivot_longer(
          cols = c(total_selected, true_positives, false_positives),
          names_to = "metric",
          values_to = "value"
        ) |>
        dplyr::mutate(
          metric = factor(
            metric,
            levels = c("total_selected", "true_positives", "false_positives"),
            labels = c("Total Selected", "True Positives", "False Positives")
          )
        )

      p <- ggplot(plot_data, aes(x = value, y = prob_factor, fill = prob_factor)) +
        ggridges::geom_density_ridges(alpha = 0.8, scale = 2) +
        facet_wrap(~metric, scales = "free_x", ncol = 3) +
        scale_fill_viridis_d(
          name = "Censoring\nprobability",
          guide = guide_legend(reverse = TRUE)
        ) +
        labs(
          title = paste0(
            "glmnet ", version, " - ",
            ifelse(alpha == 1, "Lasso", "Elastic-net"),
            " (α = ", alpha, ")"
          ),
          subtitle = "Distribution of selected variables using lambda.1se",
          x = "Number of variables",
          y = "Censoring probability"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          strip.text = element_text(size = 11, face = "bold"),
          legend.position = "right"
        )

      plot_name <- paste0(
        "v", gsub("\\.|-", "_", version), "_alpha",
        gsub("\\.", "_", as.character(alpha))
      )
      plots_list[[plot_name]] <- p
    }
  }

  plots_list
}

# Run simulation on different glmnet versions and save results ----
run_simulation <- function() {
  message("\n=== Running benchmarks with glmnet 4.1-8 ===")
  results_418 <- simulate_version("4.1-8", n_reps = 36)

  message("\n=== Running benchmarks with glmnet 4.1-9 ===")
  results_419 <- simulate_version("4.1-9", n_reps = 36)

  all_results <- rbind(results_418, results_419)
  saveRDS(all_results, "glmnet_comparison_results.rds")

  plots <- plot_densities(all_results)
  for (plot_name in names(plots)) {
    ggsave(
      filename = paste0("plot_", plot_name, ".png"),
      plot = plots[[plot_name]],
      width = 7,
      height = 7 / 1.618,
      dpi = 300
    )
  }

  all_results
}

# Line graph displaying null model rates ----
plot_null_rates <- function(results_df) {
  summary_stats <- results_df |>
    dplyr::group_by(version, alpha, prob) |>
    dplyr::summarise(
      mean_total = mean(total_selected),
      mean_tp = mean(true_positives),
      mean_fp = mean(false_positives),
      null_model_rate = mean(total_selected == 0),
      .groups = "drop"
    )

  p_null <- ggplot(
    summary_stats,
    aes(
      x = prob, y = null_model_rate,
      color = version, linetype = factor(alpha)
    )
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    scale_color_manual(values = c("4.1-8" = pal_srcery[2], "4.1-9" = pal_srcery[3])) +
    scale_linetype_manual(
      values = c("1" = "solid", "0.5" = "dashed"),
      labels = c("1" = "Lasso", "0.5" = "Elastic-net")
    ) +
    labs(
      title = "Null model rate comparison",
      subtitle = "Proportion of simulations with no variable selected",
      x = "Censoring probability",
      y = "Null model rate",
      color = "glmnet version",
      linetype = "Model type"
    ) +
    cowplot::theme_cowplot() +
    theme(
      legend.position = c(0.05, 0.70),
      legend.direction = "vertical",
      legend.key.width = unit(3, "cm")
    )

  ggsave("null_model_rate_comparison.svg", p_null, width = 7, height = 7 / 1.618)

  list(summary_stats = summary_stats, plot = p_null)
}

# Line graph comparing number of variables selected between glmnet versions ----
plot_nvar_diff <- function(results_df) {
  # Summary statistics by version and censoring probability
  summary_by_version <- results_df |>
    dplyr::group_by(version, alpha, prob) |>
    dplyr::summarise(
      mean_total_selected = mean(total_selected),
      sd_total_selected = sd(total_selected),
      median_total_selected = median(total_selected),
      null_model_rate = mean(total_selected == 0),
      mean_true_positives = mean(true_positives),
      mean_false_positives = mean(false_positives),
      sparse_model_rate = mean(total_selected <= 5),
      .groups = "drop"
    )

  difference_summary <- summary_by_version |>
    dplyr::select(alpha, prob, version, mean_total_selected, null_model_rate) |>
    tidyr::pivot_wider(names_from = version, values_from = c(mean_total_selected, null_model_rate)) |>
    dplyr::mutate(
      diff_mean_selected = `mean_total_selected_4.1-9` - `mean_total_selected_4.1-8`,
      diff_null_rate = `null_model_rate_4.1-9` - `null_model_rate_4.1-8`
    )

  # Plot difference in mean number of selected variables
  p_diff_selected <- ggplot(
    difference_summary,
    aes(
      x = prob, y = diff_mean_selected,
      color = factor(alpha), group = factor(alpha)
    )
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(
      values = c("1" = pal_srcery[1], "0.5" = pal_srcery[4]),
      labels = c("1" = "Lasso", "0.5" = "Elastic-net")
    ) +
    labs(
      title = "Change in # variables selected: glmnet 4.1-9 vs. 4.1-8",
      subtitle = "Negative values indicate fewer variables selected in glmnet 4.1-9",
      x = "Censoring probability",
      y = "Difference in mean variables selected\n(4.1-9 minus 4.1-8)",
      color = "Model type"
    ) +
    cowplot::theme_cowplot() +
    theme(
      legend.position = c(0.05, 0.85),
      legend.direction = "horizontal"
    )

  ggsave("difference_mean_selected.svg", p_diff_selected, width = 7, height = 7 / 1.618)

  list(
    summary_by_version = summary_by_version,
    difference_summary = difference_summary,
    plot = p_diff_selected
  )
}

# Heatmap for null model rates ----
plot_heatmap <- function(summary_by_version) {
  heatmap_data <- summary_by_version |>
    dplyr::mutate(
      model_type = ifelse(alpha == 1, "Lasso", "Elastic-net"),
      version_label = paste("glmnet", version)
    )

  p_heatmap <- ggplot(
    heatmap_data,
    aes(
      x = factor(prob), y = interaction(model_type, version_label),
      fill = null_model_rate
    )
  ) +
    geom_tile() +
    geom_text(
      aes(label = sprintf("%.1f%%", null_model_rate * 100)),
      color = "white", size = 3
    ) +
    scale_fill_gradient2(
      low = pal_srcery[2], mid = pal_srcery[3], high = pal_srcery[1],
      midpoint = 0.5, limits = c(0, 1),
      labels = scales::percent,
      name = "Null Model\nRate"
    ) +
    labs(
      title = "Null model rates across conditions",
      subtitle = "Proportion of simulations with no variable selection (lambda.1se)",
      x = "Censoring probability",
      y = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )

  ggsave("null_model_heatmap.svg", p_heatmap, width = 7, height = 7 / 1.618)

  p_heatmap
}

# Plot simulation outputs ----
plot_results <- function(results_df) {
  sparsity_analysis <- plot_null_rates(results_df)
  version_comparison <- plot_nvar_diff(results_df)
  heatmap_plot <- plot_heatmap(version_comparison$summary_by_version)

  invisible()
}

# Run the full analysis ----
results <- run_simulation()
plots <- plot_results(results)
