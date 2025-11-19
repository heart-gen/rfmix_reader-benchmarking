library(glue)
library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)
library(patchwork)

format_task_pop <- function(task, population_model) {
    pop_map <- c(
        two_pop = "Two Ancestries", three_pop = "Three Ancestries"
    )

    pretty_pop <- pop_map[[as.character(population_model)]]

    if (is.null(pretty_pop)) {
        pretty_pop <- population_model |>
            as.character() |> str_replace_all("_", " ") |> str_to_title()
    }

    glue("Task {task}\n({pretty_pop})")
}

make_figure <- function(df, out_prefix) {
                                        # Color-blind-friendly palette (Okabe-Ito)
    okabe_ito <- c(
        "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
    )

    df <- df |>
        filter(status == "success") |> 
        mutate(
            runtime_min = wall_time_sec / 60,
            cpu_mem_gb = peak_cpu_memory_MB / 1024,
            gpu_mem_gb = peak_gpu_memory_MB / 1024,
            task_pop = map2_chr(task, population_model, format_task_pop),
            parser_clean = factor(parser_clean, levels = sort(unique(parser_clean))),
            backend = factor(backend, levels = c("CPU", "GPU")),
            ) |>
        group_by(task_pop, backend, parser_clean) |>
        summarise(
            runtime_min_mean = mean(runtime_min, na.rm = TRUE),
            runtime_min_std  = sd(runtime_min, na.rm = TRUE),
            cpu_mem_gb_mean  = mean(cpu_mem_gb, na.rm = TRUE),
            cpu_mem_gb_std   = sd(cpu_mem_gb, na.rm = TRUE),
            gpu_mem_gb_mean  = mean(gpu_mem_gb, na.rm = TRUE),
            gpu_mem_gb_std   = sd(gpu_mem_gb, na.rm = TRUE),
            .groups = "drop"
        ) |>
        tidyr::complete(
            task_pop, backend, parser_clean,
            fill = list(
                runtime_min_mean = NA_real_,
                runtime_min_std  = NA_real_,
                cpu_mem_gb_mean  = NA_real_,
                cpu_mem_gb_std   = NA_real_,
                gpu_mem_gb_mean  = NA_real_,
                gpu_mem_gb_std   = NA_real_
            )
        )
    metric_specs <- tibble::tibble(
        metric_name = c("Runtime", "CPU memory", "GPU memory"),
        y = c("runtime_min_mean", "cpu_mem_gb_mean", "gpu_mem_gb_mean"),
        std = c("runtime_min_std", "cpu_mem_gb_std", "gpu_mem_gb_std"),
        plot_title = c("Runtime", "Peak CPU memory", "Peak GPU memory"),
        ylabel = c(
            "Runtime (min)", "Peak CPU memory (GB)", "Peak GPU memory (GB)"
        )
    )

                                        # Barplot with IQR errorbars
    plot_metric_backend <- function(data, spec, backend) {

        sub <- data |> filter(backend == !!backend)

                                        # No GPU memory for CPU runs
        if (spec$metric_name == "GPU memory" && backend == "CPU") {
            return(
                ggplot() +
                annotate("text", x = 0.5, y = 0.5,
                         label = "No GPU memory for CPU runs",
                         size = 5) +
                theme_void(base_size = 20)
            )
        }

                                        # Empty data case
        if (nrow(sub) == 0 || all(is.na(sub[[spec$y]]))) {
            return(
                ggplot() +
                annotate("text", x = 0.5, y = 0.5,
                         label = paste("No data for", backend),
                         size = 5) +
                theme_void(base_size = 20)
            )
        }

        ggplot(sub, aes(x = task_pop, y = .data[[spec$y]], fill = parser_clean)) +
            geom_col(position = position_dodge(width = 0.9)) +
            geom_errorbar(
                aes(ymin = .data[[spec$y]] - .data[[spec$std]],
                    ymax = .data[[spec$y]] + .data[[spec$std]]),
                width = 0.2, position = position_dodge(width = 0.9)
            ) +
            scale_fill_manual(values = okabe_ito) +
            labs(
                x = "", y = spec$ylabel, fill = "Parser",
                title = paste0(spec$plot_title, " - ", backend)
            ) +
            theme_bw(base_size = 20) +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "right",
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 12)
            )
    }

                                        # Build the 6 plots (3 metrics × 2 backends)
    plots <- purrr::pmap(
        metric_specs,
        \(metric_name, y, std, plot_title, ylabel) {
            spec <- list(
                metric_name = metric_name,
                y = y, std = std, plot_title = plot_title, ylabel = ylabel
            )
            list(
                cpu = plot_metric_backend(df, spec, backend = "CPU"),
                gpu = plot_metric_backend(df, spec, backend = "GPU")
            )
        }
    )

    plot_list <- unlist(plots, recursive = FALSE)

                                        # Assemble into a fixed 3x2 grid
    combined <- (plot_list[[1]] | plot_list[[2]]) /
        (plot_list[[3]] | plot_list[[4]]) /
        (plot_list[[5]] | plot_list[[6]])

    combined <- combined + patchwork::plot_layout(guides = "collect")

                                        # Save output
    png_path <- file.path(dirname(out_prefix), paste0(basename(out_prefix),
                                                      "_summary.png"))
    pdf_path <- file.path(dirname(out_prefix), paste0(basename(out_prefix),
                                                      "_summary.pdf"))

    ggsave(png_path, combined, width = 14, height = 16, dpi = 300)
    ggsave(pdf_path, combined, width = 14, height = 16)

    message("[INFO] Wrote performance figures to: ", png_path, " and ", pdf_path)

    invisible(combined)
}

make_figure_boxplot <- function(raw_df, out_prefix) {
                                        # Color-blind-friendly palette (Okabe-Ito)
    okabe_ito <- c(
        "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
    )

    df <- raw_df |>
        filter(status == "success") |> 
        mutate(
            task_pop = map2_chr(task, population_model, format_task_pop),
            parser_clean = factor(parser_clean, levels = sort(unique(parser_clean))),
            backend = factor(backend, levels = c("CPU", "GPU")),
            runtime_min = wall_time_sec / 60,
            cpu_mem_gb = peak_cpu_memory_MB / 1024,
            gpu_mem_gb = peak_gpu_memory_MB / 1024
        ) |>
        tidyr::complete(
            task_pop, backend, parser_clean,
            fill = list(
                runtime_min = NA_real_,
                cpu_mem_gb  = NA_real_,
                gpu_mem_gb  = NA_real_
            )
        )

    metric_specs <- tibble::tibble(
        metric_name = c("Runtime", "CPU memory", "GPU memory"),
        y = c("runtime_min", "cpu_mem_gb", "gpu_mem_gb"),
        plot_title = c(
            "Runtime",
            "Peak CPU memory",
            "Peak GPU memory"
        ),
        ylabel = c(
            "Runtime (min)", "Peak CPU memory (GB)", "Peak GPU memory (GB)"
        )
    )

                                        # Barplot with IQR errorbars
    plot_metric_backend <- function(data, spec, backend) {

        sub <- data |> filter(backend == !!backend)

                                        # No GPU memory for CPU runs
        if (spec$metric_name == "GPU memory" && backend == "CPU") {
            return(
                ggplot() +
                annotate("text", x = 0.5, y = 0.5,
                         label = "No GPU memory for CPU runs",
                         size = 5) +
                theme_void(base_size = 20)
            )
        }

                                        # Empty data case
        if (nrow(sub) == 0 || all(is.na(sub[[spec$y]]))) {
            return(
                ggplot() +
                annotate("text", x = 0.5, y = 0.5,
                         label = paste("No data for", backend),
                         size = 5) +
                theme_void(base_size = 20)
            )
        }
        dodge <- position_dodge(width = 0.8)
        ggplot(sub, aes(x = task_pop, y = .data[[spec$y]], fill = parser_clean)) +
            geom_jitter(position = dodge, size = 1.8, alpha = 0.6) +
            geom_boxplot(
                position = position_dodge(width = 0.8), outlier_alpha = 0.4
            ) +
            labs(
                x = "", y = spec$ylabel, fill = "Parser",
                title = paste0(spec$plot_title, " - ", backend)
            ) +
            scale_fill_manual(values = okabe_ito) +
            scale_color_manual(values = okabe_ito) + 
            theme_bw(base_size = 20) +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "right",
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 12)
            )
    }

                                        # Build the 6 plots (3 metrics × 2 backends)
    plots <- purrr::pmap(
        metric_specs,
        \(metric_name, y, plot_title, ylabel) {
            spec <- list(
                metric_name = metric_name,
                y = y,  plot_title = plot_title, ylabel = ylabel
            )
            list(
                cpu = plot_metric_backend(df, spec, backend = "CPU"),
                gpu = plot_metric_backend(df, spec, backend = "GPU")
            )
        }
    )

    plot_list <- unlist(plots, recursive = FALSE)

                                        # Assemble into a fixed 3x2 grid
    combined <- (plot_list[[1]] | plot_list[[2]]) /
        (plot_list[[3]] | plot_list[[4]]) /
        (plot_list[[5]] | plot_list[[6]])

    combined <- combined + patchwork::plot_layout(guides = "collect")

                                        # Save output
    png_path <- file.path(dirname(out_prefix), paste0(basename(out_prefix),
                                                      "_summary.png"))
    pdf_path <- file.path(dirname(out_prefix), paste0(basename(out_prefix),
                                                      "_summary.pdf"))

    ggsave(png_path, combined, width = 14, height = 16, dpi = 300)
    ggsave(pdf_path, combined, width = 14, height = 16)

    message("[INFO] Wrote performance figures to: ", png_path, " and ", pdf_path)

    invisible(combined)
}

#### Main
df <- data.table::fread("benchmark_clean.tsv")
make_figure(df, "benchmark_barplot")
make_figure_boxplot(df, "benchmark_boxplot")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
