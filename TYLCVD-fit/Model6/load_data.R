library(ggplot2)
library(dplyr)
library(rstan)
library(gridExtra)
library(outbreaks)
library(bayesplot)
library(data.table)
library(knitr)
library(kableExtra)
load_data <- function(path="/data", location="interpolated"){
    sub_path_1 <- "/home/gabrielsalcedo"
    sub_path_2 <- "Documentos/TYLCVD-fit"
    sub_path_3 <- "Model6/data"
    file_name <- "interpolated_Ramsh_data.csv"#interpolated_data#incidence_data
    data_path <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name, sep = "/")
    tomato_data <- fread(data_path, select = c("Time",
                                               "Interpolate_Ip")
                         )
    tomato_data <- data.frame(tomato_data)
#
    #reference_date <- 0#as.Date('2020-03-10')
    #final_date_sample <- 30#as.Date('2020-03-30')
    data_star_dynamics <- tomato_data# %>%
        #filter(Time >= reference_date &
        #       Time <= final_date_sample)
    head(data_star_dynamics)
    data_plot <- ggplot(data = data_star_dynamics,
                        aes(x = Time, Interpolate_Ip)) +
                    geom_bar(stat="identity", width = 0.05) +
                    geom_point() +
                    theme(axis.text.x = element_text(angle = 90)) +
                    ggtitle("LA 1582 specie data")
    #
    sub_path_1 <- "/home/gabrielsalcedo"
    sub_path_2 <- "Documentos/TYLCVD-fit"
    sub_path_3 <- "Model6/plots"
    file_name <- "tomato_input_data.pdf"
    plot_path <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name, sep = "/")
    ggsave(plot_path)
    onset <- data_star_dynamics %>%
        select(Time)
    cum_cases <- data_star_dynamics %>%
        select(Interpolate_Ip)
# 
#     cases <- data_star_dynamics %>%
#         select(Interpolated_Ip)
    # cum_cases <- unlist(cum_cases, use.names = FALSE)
    names(cum_cases)[1] <-"cum_cases"
    data <- list(onset, cum_cases, data_star_dynamics)
    return (data)
}
