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
#' plot_simulations(test)

plot_simulations <- function(adsim_df){

  # create a lollipop plot that shows mean and sd of final allele frequencies
  plot_lollipop <- function(adsim_df){
    plot_data <- adsimR::sim_overview(adsim_df)
    plot_data <- plot_data %>%
      dplyr::mutate(
        param_num = as.numeric(gsub("P", "", param_group))
      ) %>%
      dplyr::arrange(param_num) %>%
      dplyr::mutate(
        param_group = factor(param_group, levels = param_group)
      )
    P1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = c(plot_data$mean_pF, plot_data$sd_pF), y = factor(plot_data$param_group))) +
      ggplot2::geom_segment(ggplot2::aes(x = plot_data$Q1_pF, xend = plot_data$Q6_pF, y = factor(plot_data$param_group), yend = factor(plot_data$param_group), color = "95%"), linewidth = 3) +
      ggplot2::geom_segment(ggplot2::aes(x = plot_data$Q2_pF, xend = plot_data$Q5_pF, y = factor(plot_data$param_group), yend = factor(plot_data$param_group), color = "90%"), linewidth = 3) +
      ggplot2::geom_segment(ggplot2::aes(x = plot_data$Q3_pF, xend = plot_data$Q4_pF, y = factor(plot_data$param_group), yend = factor(plot_data$param_group), color = "75%"), linewidth = 3) +
      ggplot2::geom_segment(ggplot2::aes(x = plot_data$mean_pF - 0.001, xend = plot_data$mean_pF + 0.001, y = factor(plot_data$param_group), yend = factor(plot_data$param_group), color = "mean"), linewidth = 3) +
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

  # table of sim_overview output saved as grob
  p_table <- function(adsim_df){
    table_data <- adsimR::sim_overview(adsim_df) %>%
      dplyr::select(-adsim_df$Q1_pF, -adsim_df$Q2_pF, -adsim_df$Q3_pF, -adsim_df$Q4_pF, -adsim_df$Q5_pF, -adsim_df$Q6_pF)
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
  if (nrow(table4grob) < 15) {
    return(gridExtra::grid.arrange(grobs = list(plot4grob, table4grob), ncol = 1, width = c(2,1)))
  } else {
    return(plot4grob)
  }

}
