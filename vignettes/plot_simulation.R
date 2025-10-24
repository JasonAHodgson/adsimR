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
  plot_data <- overview_adsim(test)

 P1 <- ggplot(plot_data, aes(x = c(mean_pF, sd_pF), y = factor(param_group))) +
    #geom_segment(aes(x = mean_pF - sd_pF, xend=mean_pF + sd_pF, y=factor(param_group), yend=factor(param_group)), color = "blue", linewidth=13) +
    geom_segment(aes(x = Q1_pF, xend=Q6_pF, y=factor(param_group), yend=factor(param_group)), color = "#AA336A", linewidth=13) +
    geom_segment(aes(x = Q2_pF, xend=Q5_pF, y=factor(param_group), yend=factor(param_group)), color = "#CF9FFF", linewidth=13) +
    geom_segment(aes(x = Q3_pF, xend=Q4_pF, y=factor(param_group), yend=factor(param_group)), color = "#E6E6FA", linewidth=13) +
    geom_segment(aes(x = mean_pF - 0.001, xend=mean_pF + 0.001, y=factor(param_group), yend=factor(param_group)), color = "#702963", linewidth=13) +
    labs(title = "Allele Frequency Summary (Mean ± SD)",
         x = "Allele Frequency",
         y = "Parameter Set"
    )+
    theme_minimal()
 return(P1)
}
plot_lollipop(test)
