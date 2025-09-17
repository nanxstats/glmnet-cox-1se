# Compare lambda.1se rule behavior for Cox regressions in glmnet 4.1-8 vs. 4.1-9
# Examine the effect of CV error normalization changes on model selection

# Use color palette
pal_obs <- ggsci::pal_observable()(6)

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

  `%dofuture%` <- doFuture::`%dofuture%`

  future::plan(future::multisession, workers = parallelly::availableCores() - 1)

  set.seed(42)
  results_list <- foreach::foreach(
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

  future::plan(future::sequential)

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

  # Create plot for each alpha (model type) with both versions
  plots_list <- list()

  for (alpha in unique(results_df$alpha)) {
    # Prepare data with version info for faceting
    plot_data <- results_df |>
      dplyr::filter(alpha == !!alpha) |>
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
        ),
        version_label = paste0("glmnet ", version)
      )

    # Define x-axis limits for each metric and model type
    x_limits <- list()
    if (alpha == 0.5) {
      x_limits[["Total Selected"]] <- c(-1, 40)
      x_limits[["True Positives"]] <- c(-1, 13)
      x_limits[["False Positives"]] <- c(-1, 20)
    }
    if (alpha == 1) {
      x_limits[["Total Selected"]] <- c(-1, 30)
      x_limits[["True Positives"]] <- c(-1, 13)
      x_limits[["False Positives"]] <- c(-1, 10)
    }

    # Create individual panels for cowplot alignment
    panel_list <- list()

    for (ver in unique(plot_data$version)) {
      for (met in levels(plot_data$metric)) {
        panel_data <- plot_data |>
          dplyr::filter(version == ver, metric == met)

        p_panel <- ggplot2::ggplot(panel_data, ggplot2::aes(x = value, y = prob_factor, fill = prob_factor)) +
          ggridges::geom_density_ridges(alpha = 0.8, scale = 2) +
          ggplot2::scale_x_continuous(limits = x_limits[[met]]) +
          ggplot2::scale_fill_viridis_d(
            name = "Censoring\nprobability",
            guide = ggplot2::guide_legend(reverse = TRUE),
            option = "plasma",
            direction = -1
          ) +
          ggplot2::labs(
            title = paste0(panel_data$version_label[1], " - ", met),
            x = "Number of variables",
            y = if (met == "Total Selected") "Censoring probability" else ""
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
            axis.title.y = ggplot2::element_text(size = 10),
            axis.title.x = ggplot2::element_text(size = 10),
            legend.position = "none"
          )

        panel_name <- paste0(ver, "_", gsub(" ", "_", met))
        panel_list[[panel_name]] <- p_panel
      }
    }

    # Extract legend from one plot
    legend_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value, y = prob_factor, fill = prob_factor)) +
      ggridges::geom_density_ridges(alpha = 0.8, scale = 2) +
      ggplot2::scale_fill_viridis_d(
        name = "Censoring\nprobability",
        guide = ggplot2::guide_legend(reverse = TRUE),
        option = "plasma",
        direction = -1
      ) +
      ggplot2::theme_minimal()
    legend <- cowplot::get_legend(legend_plot)

    # Arrange panels with aligned axes using cowplot
    aligned_panels <- cowplot::align_plots(
      panel_list[["4.1-8_Total_Selected"]], panel_list[["4.1-8_True_Positives"]], panel_list[["4.1-8_False_Positives"]],
      panel_list[["4.1-9_Total_Selected"]], panel_list[["4.1-9_True_Positives"]], panel_list[["4.1-9_False_Positives"]],
      align = "hv", axis = "tblr"
    )

    # Create the combined plot
    combined_plot <- cowplot::plot_grid(
      aligned_panels[[1]], aligned_panels[[2]], aligned_panels[[3]],
      aligned_panels[[4]], aligned_panels[[5]], aligned_panels[[6]],
      ncol = 3, nrow = 2,
      rel_widths = c(1, 1, 1),
      rel_heights = c(1, 1)
    )

    # Add overall title and legend
    title <- cowplot::ggdraw() +
      cowplot::draw_label(
        paste0(
          ifelse(alpha == 1, "Lasso", "Elastic-net"),
          " (α = ", alpha, ") - Distribution of selected variables using lambda.1se"
        ),
        fontface = "bold",
        size = 14,
        x = 0.5,
        hjust = 0.5
      )

    # Combine title, plot, and legend
    final_plot <- cowplot::plot_grid(
      title,
      cowplot::plot_grid(combined_plot, legend, ncol = 2, rel_widths = c(6, 1)),
      ncol = 1,
      rel_heights = c(0.1, 1)
    )

    plot_name <- paste0("alpha", gsub("\\.", "-", as.character(alpha)))
    plots_list[[plot_name]] <- final_plot
  }

  plots_list
}

# Run simulation on different glmnet versions and save results ----
run_simulation <- function() {
  message("\n=== Running benchmarks with glmnet 4.1-8 ===")
  results_418 <- simulate_version("4.1-8", n_reps = 500)

  message("\n=== Running benchmarks with glmnet 4.1-9 ===")
  results_419 <- simulate_version("4.1-9", n_reps = 500)

  all_results <- rbind(results_418, results_419)
  saveRDS(all_results, "glmnet_comparison_results.rds")

  plots <- plot_densities(all_results)
  for (plot_name in names(plots)) {
    ggplot2::ggsave(
      filename = paste0("plot-", plot_name, ".png"),
      plot = plots[[plot_name]],
      width = 12,
      height = 8,
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

  p_null <- ggplot2::ggplot(
    summary_stats,
    ggplot2::aes(
      x = prob, y = null_model_rate,
      color = version, linetype = factor(alpha)
    )
  ) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = c("4.1-8" = pal_obs[1], "4.1-9" = pal_obs[4])) +
    ggplot2::scale_linetype_manual(
      values = c("1" = "solid", "0.5" = "dashed"),
      labels = c("1" = "Lasso", "0.5" = "Elastic-net")
    ) +
    ggplot2::labs(
      title = "Null model rate comparison",
      subtitle = "Proportion of simulations with no variable selected",
      x = "Censoring probability",
      y = "Null model rate",
      color = "glmnet version",
      linetype = "Model type"
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.05, 0.70),
      legend.direction = "vertical",
      legend.key.width = grid::unit(1.3, "cm")
    )

  ggplot2::ggsave("null-model-rate-comparison.svg", p_null, width = 7, height = 7 / 1.618)

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
  p_diff_selected <- ggplot2::ggplot(
    difference_summary,
    ggplot2::aes(x = prob, y = diff_mean_selected, color = factor(alpha), group = factor(alpha))
  ) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::scale_color_manual(
      values = c("1" = pal_obs[3], "0.5" = pal_obs[5]),
      labels = c("1" = "Lasso", "0.5" = "Elastic-net")
    ) +
    ggplot2::labs(
      title = "Change in # variables selected: glmnet 4.1-9 vs. 4.1-8",
      subtitle = "Negative values indicate fewer variables selected in glmnet 4.1-9",
      x = "Censoring probability",
      y = "Difference in mean variables selected\n(4.1-9 minus 4.1-8)",
      color = "Model type"
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.05, 0.85),
      legend.direction = "horizontal",
      legend.key.width = grid::unit(1.3, "cm")
    )

  ggplot2::ggsave("difference-mean-selected.svg", p_diff_selected, width = 7, height = 7 / 1.618)

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

  p_heatmap <- ggplot2::ggplot(
    heatmap_data,
    ggplot2::aes(
      x = factor(prob), y = interaction(model_type, version_label),
      fill = null_model_rate
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", null_model_rate * 100)),
      color = "white", size = 3
    ) +
    ggplot2::scale_fill_gradient2(
      low = pal_obs[5], mid = pal_obs[2], high = pal_obs[3],
      midpoint = 0.5, limits = c(0, 1),
      labels = scales::percent,
      name = "Null Model\nRate"
    ) +
    ggplot2::labs(
      title = "Null model rates across conditions",
      subtitle = "Proportion of simulations with no variable selection (lambda.1se)",
      x = "Censoring probability",
      y = ""
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12),
      axis.text.y = ggplot2::element_text(size = 10),
      legend.position = "right"
    )

  ggplot2::ggsave("null-model-heatmap.svg", p_heatmap, width = 7, height = 7 / 1.618)

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
