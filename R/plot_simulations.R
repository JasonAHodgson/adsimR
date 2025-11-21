#' plot_simulations
#'
#' This function creates an autoplot of adsim output
#'
#' @param adsim_df  A data frame containing three columns:
#' - Sim: a numeric index containing the simulation run
#' - p0: a numeric containing the starting allele following the first
#'   generation of admixture.
#' - pF: a numeric containing the final allele frequency of each simulation.
#' @return a plot containing two panes:
#' - a plot containing the summary of allele frequencies of the simulation,
#' including quartiles at 75%, 90%, and 95%.
#' - a table containing the parameters of each simulation.
#' @export
#' @examples
#' test <- adsim(ngens = 10, k = 2, l = 5, admix = c(0.2, 0.5),
#' m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
#' p2 = 0.9, nsims = 10)
#'
#' plot_func(test)

plot_simulations <- function(adsim_df){
  overview_adsim <- function(adsim_df){
    sum_adsim_df <- adsim_df %>%
      dplyr::group_by(ngens, k, l, admix, m1, m2, ne, p1, p2) %>%
      dplyr::summarise(
        mean_pF = mean(pF),
        sd_pF   = stats::sd(pF),
        # calculate percentiles
        Q1_pF = stats::quantile(pF, 0.05),
        Q2_pF = stats::quantile(pF, 0.10),
        Q3_pF = stats::quantile(pF, 0.25),
        Q4_pF = stats::quantile(pF, 0.75),
        Q5_pF = stats::quantile(pF, 0.90),
        Q6_pF = stats::quantile(pF, 0.95),
        .groups = "drop"
      )
    # Create vector P1..PN where N = number of rows
    sum_adsim_df$param_group <- paste0("P", seq_len(nrow(sum_adsim_df)))
    return(sum_adsim_df)
  }

  # create a lollipop plot that shows mean and sd of final allele frequencies
  plot_lollipop <- function(adsim_df){
    plot_data <- overview_adsim(adsim_df)

    P1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = c(mean_pF, sd_pF), y = factor(param_group))) +
      ggplot2::geom_segment(ggplot2::aes(x = Q1_pF, xend = Q6_pF, y = factor(param_group), yend = factor(param_group), color = "95%"), linewidth = 13) +
      ggplot2::geom_segment(ggplot2::aes(x = Q2_pF, xend = Q5_pF, y = factor(param_group), yend = factor(param_group), color = "90%"), linewidth = 13) +
      ggplot2::geom_segment(ggplot2::aes(x = Q3_pF, xend = Q4_pF, y = factor(param_group), yend = factor(param_group), color = "75%"), linewidth = 13) +
      ggplot2::geom_segment(ggplot2::aes(x = mean_pF - 0.001, xend = mean_pF + 0.001, y = factor(param_group), yend = factor(param_group), color = "mean"), linewidth = 13) +
      ggplot2::scale_color_manual(name = "Legend",
                                  values = c("95%" = "#AA336A",
                                             "90%" = "#CF9FFF",
                                             "75%" = "#E6E6FA",
                                             "mean" = "#702963")) +
      ggplot2::labs(title = "Allele Frequency Summary",
                    x = "Allele Frequency",
                    y = "Parameter Set") +
      ggplot2::theme_minimal()
    return(P1)
  }
  plot4grob <- plot_lollipop(adsim_df)

  # table of overview_adsim output saved as grob
  p_table <- function(adsim_df){
    table_data <- overview_adsim(adsim_df) %>%
      dplyr::select(-Q1_pF, -Q2_pF, -Q3_pF, -Q4_pF, -Q5_pF, -Q6_pF)
    theme <- gridExtra::ttheme_default(
      core = list(
        bg_params = list(fill = rep("white", nrow(adsim_df)))
      ),
      colhead = list(
        fg_params = list(col = "black", fontface = "bold"),
        bg_params = list(fill = rep("#c9d", ncol(table_data)))
      )
    )
    table_grob <- gridExtra::tableGrob(table_data, theme = theme)
    return(table_grob)
  }
  table4grob <- p_table(adsim_df)

  # combine plot and table into single grob
  return(gridExtra::grid.arrange(plot4grob, table4grob, ncol = 1))
}

