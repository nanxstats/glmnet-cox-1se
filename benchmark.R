# Compare Cox regression lambda.1se rule behavior in glmnet 4.1-8 vs. 4.1-9
# Examine the effect of CV error normalization changes on model selection

library(glmnet)
library(doFuture)
library(ggplot2)

# Simulate Cox regression ----
simulate_cox <- function(
    n = 1000, p = 100, nzc = 10,
    prob, alpha,
    use_1se = TRUE) {
  # Generate data
  x <- matrix(rnorm(n * p), n, p)

  # True coefficients
  beta <- rnorm(nzc)
  fx <- x[, seq(nzc)] %*% beta / 3
  hx <- exp(fx)
  ty <- rexp(n, hx)

  # Censoring (prob is the probability of being censored)
  tcens <- rbinom(n = n, prob = prob, size = 1)
  y <- cbind(time = ty, status = 1 - tcens)

  # Run cross-validation
  fit_cv <- cv.glmnet(x, y, family = "cox", nfolds = 5, alpha = alpha)

  # Select lambda
  lambda_selected <- if (use_1se) fit_cv$lambda.1se else fit_cv$lambda.min

  # Fit model on complete data with selected lambda
  fit <- glmnet(x, y, family = "cox", alpha = alpha, lambda = lambda_selected)

  # Calculate metrics
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

# Simulate for a specific glmnet version ----
simulate_version <- function(version, n_reps = 500) {
  # Set up parameters
  prob_values <- seq(0.1, 0.9, by = 0.1)
  alpha_values <- c(1, 0.5) # lasso and elastic-net

  # Forcefully unload glmnet if loaded
  if ("glmnet" %in% loadedNamespaces()) unloadNamespace("glmnet")

  # Install specific version (force ensures we get the exact version)
  message(paste("Installing glmnet", version, "..."))
  remotes::install_version("glmnet", version, quiet = TRUE, upgrade = "never", force = TRUE)

  # Load the specific version
  library(glmnet)

  # Verify version
  current_version <- packageVersion("glmnet")
  message(paste("Running experiments with glmnet version:", current_version))

  # Set up parallel backend
  plan(multisession, workers = parallelly::availableCores() - 1)

  # Set seed for reproducibility
  set.seed(42)

  # Run experiments using lambda.1se (main focus)
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

  # Reset to sequential processing
  plan(sequential)

  # Combine results using dplyr
  results <- dplyr::bind_rows(results_list)

  as.data.frame(results)
}

# Plot selection densities ----
plot_densities <- function(results_df) {
  # Convert prob to factor for proper ordering
  results_df$prob_factor <- factor(
    results_df$prob,
    levels = seq(0.1, 0.9, by = 0.1),
    ordered = TRUE
  )

  # Create plot for each version and alpha combination
  plots_list <- list()

  for (version in unique(results_df$version)) {
    for (alpha in unique(results_df$alpha)) {
      # Filter data
      plot_data <- results_df |>
        dplyr::filter(version == !!version, alpha == !!alpha) |>
        tidyr::pivot_longer(
          cols = c(total_selected, true_positives, false_positives),
          names_to = "metric",
          values_to = "value"
        ) |>
        dplyr::mutate(
          metric = factor(metric,
            levels = c("total_selected", "true_positives", "false_positives"),
            labels = c("Total Selected", "True Positives", "False Positives")
          )
        )

      # Create ridgeline plot
      p <- ggplot(plot_data, aes(x = value, y = prob_factor, fill = prob_factor)) +
        ggridges::geom_density_ridges(alpha = 0.8, scale = 2) +
        facet_wrap(~metric, scales = "free_x", ncol = 3) +
        scale_fill_viridis_d(name = "Censoring\nProbability", direction = -1) +
        labs(
          title = paste0(
            "glmnet ", version, " - ",
            ifelse(alpha == 1, "Lasso", "Elastic-Net"),
            " (α = ", alpha, ")"
          ),
          subtitle = "Distribution of selected variables using lambda.1se",
          x = "Number of Variables",
          y = "Censoring Probability"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          strip.text = element_text(size = 11, face = "bold"),
          legend.position = "right"
        )

      # Store plot
      plot_name <- paste0(
        "v", gsub("\\.", "_", version), "_alpha",
        gsub("\\.", "_", as.character(alpha))
      )
      plots_list[[plot_name]] <- p
    }
  }

  plots_list
}

# Compare glmnet versions ----
compare_glmnet_versions <- function() {
  message("Starting glmnet version comparison experiments...")

  # Run experiments for both versions
  message("\n=== Running experiments for glmnet 4.1-8 ===")
  results_418 <- simulate_version("4.1-8", n_reps = 27)

  message("\n=== Running experiments for glmnet 4.1-9 ===")
  results_419 <- simulate_version("4.1-9", n_reps = 27)

  # Combine results
  all_results <- rbind(results_418, results_419)

  # Save results
  saveRDS(all_results, "glmnet_comparison_results.rds")
  message("\nResults saved to glmnet_comparison_results.rds")

  # Create plots
  message("\nCreating density visualizations...")
  plots <- plot_densities(all_results)

  # Save plots
  for (plot_name in names(plots)) {
    ggsave(
      filename = paste0("plot_", plot_name, ".png"),
      plot = plots[[plot_name]],
      width = 12,
      height = 8,
      dpi = 300
    )
  }

  message("\nPlots saved as PNG files")

  all_results
}

# Summarize sparsity ----
summarize_sparsity <- function(results_df) {
  # Calculate summary statistics
  summary_stats <- results_df |>
    dplyr::group_by(version, alpha, prob) |>
    dplyr::summarise(
      mean_total = mean(total_selected),
      mean_tp = mean(true_positives),
      mean_fp = mean(false_positives),
      null_model_rate = mean(total_selected == 0),
      .groups = "drop"
    )

  # Plot null model rates
  p_null <- ggplot(
    summary_stats,
    aes(
      x = prob, y = null_model_rate,
      color = version, linetype = factor(alpha)
    )
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    scale_color_manual(values = c("4.1-8" = "#E69F00", "4.1-9" = "#56B4E9")) +
    scale_linetype_manual(
      values = c("1" = "solid", "0.5" = "dashed"),
      labels = c("1" = "Lasso", "0.5" = "Elastic-Net")
    ) +
    labs(
      title = "Null Model Rate Comparison",
      subtitle = "Proportion of simulations resulting in no variable selection",
      x = "Censoring Probability",
      y = "Null Model Rate",
      color = "glmnet Version",
      linetype = "Model Type"
    ) +
    cowplot::theme_cowplot() +
    theme(legend.position = "bottom")

  ggsave("null_model_rate_comparison.png", p_null, width = 10, height = 6, dpi = 300)

  list(summary_stats = summary_stats, plot = p_null)
}

# Compare versions ----
compare_versions <- function(results_df) {
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
      values = c("1" = "#E69F00", "0.5" = "#56B4E9"),
      labels = c("1" = "Lasso", "0.5" = "Elastic-Net")
    ) +
    labs(
      title = "Change in Variable Selection: glmnet 4.1-9 vs 4.1-8",
      subtitle = "Negative values indicate fewer variables selected in version 4.1-9",
      x = "Censoring Probability",
      y = "Difference in Mean Variables Selected\n(4.1-9 minus 4.1-8)",
      color = "Model Type"
    ) +
    cowplot::theme_cowplot() +
    theme(legend.position = "bottom")

  ggsave("difference_mean_selected.png", p_diff_selected, width = 10, height = 6, dpi = 300)

  list(
    summary_by_version = summary_by_version,
    difference_summary = difference_summary,
    plot = p_diff_selected
  )
}

# Plot heatmap ----
plot_heatmap <- function(summary_by_version) {
  heatmap_data <- summary_by_version |>
    dplyr::mutate(
      model_type = ifelse(alpha == 1, "Lasso", "Elastic-Net"),
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
      low = "darkgreen", mid = "yellow", high = "darkred",
      midpoint = 0.5, limits = c(0, 1),
      labels = scales::percent,
      name = "Null Model\nRate"
    ) +
    labs(
      title = "Null Model Rates Across Conditions",
      subtitle = "Proportion of simulations with no variable selection (lambda.1se)",
      x = "Censoring Probability",
      y = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )

  ggsave("null_model_heatmap.png", p_heatmap, width = 12, height = 6, dpi = 300)

  p_heatmap
}

# Analyze results ----
analyze_results <- function(results_df) {
  sparsity_analysis <- summarize_sparsity(results_df)
  version_comparison <- compare_versions(results_df)
  heatmap_plot <- plot_heatmap(version_comparison$summary_by_version)

  invisible()
}

# Run the full analysis pipeline
results <- compare_glmnet_versions()
analysis_results <- analyze_results(results)
