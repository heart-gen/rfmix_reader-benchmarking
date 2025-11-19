library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(glue)
library(patchwork)

format_task_pop <- function(task, population_model) {
  pop_map <- c(
    two_pop = "Two Ancestries",
    three_pop = "Three Ancestries"
  )

  pretty_pop <- pop_map[[as.character(population_model)]]

  if (is.null(pretty_pop)) {
    pretty_pop <- population_model |>
      as.character() |>
      str_replace_all("_", " ") |>
      str_to_title()
  }

  glue("Task {task}\n({pretty_pop})")
}

make_figure_fixed <- function(summary, out_prefix) {

  df <- summary |>
    mutate(
      task_pop = map2_chr(task, population_model, format_task_pop),
      backend = factor(backend, levels = c("CPU", "GPU"))
    )

  metric_specs <- tibble::tibble(
    metric_name = c("Runtime", "CPU memory", "GPU memory"),
    y = c(
      "wall_time_median_sec",
      "cpu_mem_median_mb",
      "gpu_mem_median_mb"
    ),
    q1 = c(
      "wall_time_q1_sec",
      "cpu_mem_q1_mb",
      "gpu_mem_q1_mb"
    ),
    q3 = c(
      "wall_time_q3_sec",
      "cpu_mem_q3_mb",
      "gpu_mem_q3_mb"
    ),
    plot_title = c(
      "A. Runtime (median [IQR])",
      "B. Peak CPU memory (median [IQR])",
      "C. Peak GPU memory (median [IQR])"
    ),
    ylabel = c(
      "Runtime (s)",
      "Peak CPU memory (MB)",
      "Peak GPU memory (MB)"
    )
  )

  # ---- Utility: barplot with IQR errorbars ----
  plot_metric_backend <- function(data, spec, backend) {

    sub <- data |> filter(backend == !!backend)

    # No GPU memory for CPU runs
    if (spec$metric_name == "GPU memory" && backend == "CPU") {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = "No GPU memory for CPU runs",
                   size = 5) +
          theme_void(base_size = 15)
      )
    }

    # Empty data case
    if (nrow(sub) == 0 || all(is.na(sub[[spec$y]]))) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = paste("No data for", backend),
                   size = 5) +
          theme_void(base_size = 15)
      )
    }

    ggplot(sub, aes(x = task_pop, y = .data[[spec$y]], fill = parser_clean)) +
      geom_col(position = position_dodge(width = 0.9)) +
      geom_errorbar(
        aes(
          ymin = .data[[spec$y]] - (.data[[spec$y]] - .data[[spec$q1]]),
          ymax = .data[[spec$y]] + (.data[[spec$q3]] - .data[[spec$y]])
        ),
        width = 0.2,
        position = position_dodge(width = 0.9)
      ) +
      labs(
        x = "",
        y = spec$ylabel,
        title = paste0(spec$plot_title, " - ", backend),
        fill = "Parser"
      ) +
      theme_bw(base_size = 15) +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = if (backend == "GPU") "right" else "none"
      )
  }

  # ---- Build the 6 plots (3 metrics × 2 backends) ----
  plots <- purrr::pmap(
    metric_specs,
    \(metric_name, y, q1, q3, plot_title, ylabel) {
      spec <- list(
        metric_name = metric_name,
        y = y, q1 = q1, q3 = q3,
        plot_title = plot_title,
        ylabel = ylabel
      )
      list(
        cpu = plot_metric_backend(df, spec, backend = "CPU"),
        gpu = plot_metric_backend(df, spec, backend = "GPU")
      )
    }
  )

  # Flatten to a list of length 6 in row-major order (Runtime row, CPU→GPU)
  plot_list <- unlist(plots, recursive = FALSE)

  # Assemble into a fixed 3x2 grid
  combined <- (plot_list[[1]] | plot_list[[2]]) /
              (plot_list[[3]] | plot_list[[4]]) /
              (plot_list[[5]] | plot_list[[6]])

  # ---- Save output ----
  png_path <- file.path(dirname(out_prefix), paste0(basename(out_prefix), "_summary.png"))
  pdf_path <- file.path(dirname(out_prefix), paste0(basename(out_prefix), "_summary.pdf"))

  ggsave(png_path, combined, width = 14, height = 16, dpi = 300)
  ggsave(pdf_path, combined, width = 14, height = 16)

  message("[INFO] Wrote performance figures to: ", png_path, " and ", pdf_path)

  invisible(combined)
}
