library(dplyr)
library(tidyr)
test <- adsim(ngens = 10, k = 2, l = 5, admix = c(0.2, 0.5),
              m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
              p2 = 0.9, nsims = 10)
overview_adsim <- function(df){
  sum_df <- df %>%
    group_by(ngens, k, l, admix, m1, m2, ne, p1, p2) %>%
    summarise(
      mean_pF = mean(pF),
      sd_pF = sd(pF),
    #calculate percentiles
      Q1_pF = quantile(pF, 0.05),
      Q2_pF = quantile(pF, 0.10),
      Q3_pF = quantile(pF, 0.25),
      Q4_pF = quantile(pF, 0.75),
      Q5_pF = quantile(pF, 0.90),
      Q6_pF = quantile(pF, 0.95),
      .groups = "drop"
    )
    # Create vector P1..PN where N = number of rows
    sum_df$param_group <- paste0("P", seq_len(nrow(sum_df)))
  return(sum_df)
}
overview_adsim(test)


library(ggplot2)
#create a lollipop plot that would show the mean and sd of the final allele frequencies for each parameter set
plot_lollipop <- function(df){
  plot_data <- overview_adsim(df)

 P1 <- ggplot(plot_data, aes(x = c(mean_pF, sd_pF), y = factor(param_group))) +
    geom_segment(aes(x = Q1_pF, xend=Q6_pF, y=factor(param_group), yend=factor(param_group), color = "95%"), linewidth=13) +
    geom_segment(aes(x = Q2_pF, xend=Q5_pF, y=factor(param_group), yend=factor(param_group), color = "90%"), linewidth=13) +
    geom_segment(aes(x = Q3_pF, xend=Q4_pF, y=factor(param_group), yend=factor(param_group), color = "75%"), linewidth=13) +
    geom_segment(aes(x = mean_pF - 0.001, xend=mean_pF + 0.001, y=factor(param_group), yend=factor(param_group), color = "mean"), linewidth=13) +
    scale_color_manual(name = "Legend",
                       values = c("95%" = "#AA336A",
                                  "90%" = "#CF9FFF",
                                  "75%" = "#E6E6FA",
                                  "mean" = "#702963"))+
    labs(title = "Allele Frequency Summary (Mean ± SD)",
         x = "Allele Frequency",
         y = "Parameter Set"
    )+
    theme_minimal()
 return(P1)
}
plot4grob <- plot_lollipop(test)


#A table of the overview_adsim output and saved as a grob
library(knitr)
library(kableExtra)
p_table <- function(df){
  table_data <- overview_adsim(df) %>%
    select(-Q1_pF, -Q2_pF, -Q3_pF, -Q4_pF, -Q5_pF, -Q6_pF)
  theme <- ttheme_default(
    core = list(
      bg_params = list(fill = c(rep("white", nrow(df)))
    )),
    colhead = list(
      fg_params = list(col = "black", fontface = "bold"),
      bg_params = list(fill = c(rep("#c9d", ncol(table_data))))
      )
  )
  table_grob <- tableGrob(table_data, theme = theme)
  return(table_grob)
}
table4grob <- p_table(test)


#Create a single grob that contains both the plot and the table
library(gridExtra)
grid.arrange(plot4grob, table4grob, ncol = 1)

