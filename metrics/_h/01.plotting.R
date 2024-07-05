## This script plots the metrics comparison

library(dplyr)
library(ggpubr)

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

load_data <- function(COMPUTE){
    df <- data.table::fread("../_h/extracted_metrics.csv") |>
        mutate(name=ifelse(make_binaries,
                           "RFMix-reader\n(Create Binaries)", method)) |>
        select(-method, -make_binaries)
    if(COMPUTE){
        df <- select(df, -mem)
    } else {
        df <- select(df, -compute_mem)
    }
    return(tidyr::pivot_longer(df, !c(name, processor, pop, simu),
                               names_to="metrics", values_to="value") |>
           data.table::as.data.table() |>
           mutate_if(is.character, as.factor))
}

plotting <- function(df){
    barplt <- ggbarplot(df, x="metrics", y="value", fill="name",
                        position=position_dodge(0.9), label=TRUE,
                        palette="npg", facet.by="processor")
    save_ggplots("testing", barplt, 12, 7)
}

#### Main

real_df <- filter(load_data(TRUE), !simu) |>
            select(name, processor, metrics, value)
simu_df <- filter(load_data(), simu) |>
            select(name, pop, processor, metrics, value)
