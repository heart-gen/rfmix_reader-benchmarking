## This script plots the metrics comparison

library(dplyr)
library(ggpubr)

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

load_data <- function(){
    df <- data.table::fread("../_h/extracted_metrics.csv") |>
           mutate(name=ifelse(make_binaries,
                              "RFMix-reader\n(Create Binaries)",
                              method),
                  time = round(time / 60, 1),
                  mem = round(mem / 1000, 1),
                  compute_mem = round(compute_mem / 1000, 1)) |>
           select(-method, -make_binaries) |>
        mutate_if(is.character, as.factor)
    return(df)
}

plotting_simu <- function(metric, PROCESSOR, ylab){
    outfile <- paste(tolower(PROCESSOR), "barplot_simulation",
                     metric, sep=".")
    df <- filter(load_data(), simu) |>
        filter(processor == PROCESSOR)
    if(PROCESSOR == "CPU"){
        df$name <- factor(df$name, levels=c("pandas",
                                            "RFMix-reader\n(Create Binaries)",
                                            "RFMix-reader"))
    } else {
        df$name <- factor(df$name, levels=c("cuDF",
                                            "RFMix-reader\n(Create Binaries)",
                                            "RFMix-reader"))
    }
    barplt <- ggbarplot(df, x="name", y=metric, fill="name",
                        width=0.5, label=TRUE, palette="grey",
                        facet.by="pop", ylab=ylab, xlab="",
                        ggtheme=theme_pubr(base_size=12, border=TRUE)) +
        rotate_x_text(45) + font("ylab", size=15) + rremove("legend")
    save_ggplots(outfile, barplt, 5, 6)
}

plotting_real <- function(metric, PROCESSOR, ylab){
    outfile <- paste(tolower(PROCESSOR), "barplot_real",
                     metric, sep=".")
    df <- filter(load_data(), !simu) |>
        filter(processor == PROCESSOR)
    if(PROCESSOR == "CPU"){
        df$name <- factor(df$name, levels=c("pandas",
                                            "RFMix-reader\n(Create Binaries)",
                                            "RFMix-reader"))
    } else {
        df$name <- factor(df$name, levels=c("cuDF",
                                            "RFMix-reader\n(Create Binaries)",
                                            "RFMix-reader"))
    }
    barplt <- ggbarplot(df, x="name", y=metric, fill="name", width=0.5,
                        label=TRUE, palette="grey", ylab=ylab, xlab="",
                        ggtheme=theme_pubr(base_size=12, border=TRUE)) +
        rotate_x_text(45) + font("ylab", size=15) + rremove("legend")
    save_ggplots(outfile, barplt, 3, 6)
}

#### Main

for(processor in c("CPU", "GPU")){
    plotting_simu("time", processor, "Computational Time (minutes)")
    plotting_real("time", processor, "Computational Time (minutes)")
    plotting_simu("mem", processor, "Memory (GB)")
    plotting_real("mem", processor, "Memory (GB)")
    plotting_simu("compute_mem", processor, "Compute Memory (GB)")
    plotting_real("compute_mem", processor, "Compute Memory (GB)")
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
