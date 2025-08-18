#' adsimR
#'
#' This function simulates evolution of an allele through genetic drift in an
#' admixed population. The user must specify the starting allele frequency in
#' both parent populations, the number of generations to simulate, the starting
#' effective population size, the admixture proportion, the population growth
#' parameters, the migration rate from both parent populations, and the number
#' of simulations to perform. The function returns a data frame containing the
#' final allele frequency of each simulation.
#'
#' @param ngens Number of generations to simulate.
#' @param k population growth parameter.
#' @param l number of generations prior to population growth.
#' @param admix Proportion of ancestry contributed by the first parent population.
#' @param ne effective population size of generation 0.
#' @param p1 Initial allele frequency in the first parent population.
#' @param p2 Initial allele frequency in the second parent population.
#' @param m1 Migration rate from the first parent population.
#' @param m2 Migration rate from the second parent population.
#' @param nsims Number of simulations to perform.


