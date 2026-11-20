#' Populations of the 30 largest cities of the United States in 1940
#'
#' The example of Hill (1975, Table 1), who took the populations from the
#' Statistical Abstract of the United States, 1950, p. 57; cities are defined
#' by their political boundaries. The table prints the populations in
#' thousands; they are stored here in inhabitants, the unit in which the
#' paper's estimates of `k C` are reproduced.
#'
#' @format A data frame with 30 rows: `rank` and `population` (inhabitants,
#'   rounded to the thousand or, for ranks 28-30, to the hundred).
#' @source Hill, B. M. (1975). A simple general approach to inference about
#'   the tail of a distribution. The Annals of Statistics, 3(5), 1163-1174,
#'   Table 1.
#' @examples
#' hill1975(us_cities_1940$population, k = 29)
"us_cities_1940"

#' Exceedances of the tensile strength of nylon carpet fibres
#'
#' The example of Grimshaw (1993, Table 1): the amounts (kg/mm^2) by which
#' the tensile strength of 15 nylon carpet fibres exceeded a testing
#' threshold, which is not given because it is proprietary. The note of the
#' table prints a median of .16200, but the printed values have median .184;
#' the other summaries of the note (mean .25267, standard deviation .24360,
#' interquartile range .29000) agree with the data.
#'
#' @format A data frame with 15 rows and one column, `excess`, in the order
#'   of the table.
#' @source Grimshaw, S. D. (1993). Computing maximum likelihood estimates for
#'   the generalized Pareto distribution. Technometrics, 35(2), 185-191,
#'   Table 1.
#' @examples
#' grimshaw1993(nylon_fibres$excess)
"nylon_fibres"

#' Areas of the largest Australian bushfires of 1986-1987
#'
#' The areas of the 30 largest of 75 bushfires recorded in Australia in
#' 1986-1987, the first example of Wang (1995), who took them from Smith
#' (1993).
#'
#' @format A data frame with 30 rows and one column, `area` (hectares), in
#'   decreasing order.
#' @source Wang, J. Z. (1995). Selection of the k largest order statistics
#'   for the domain of attraction of the Gumbel distribution. Journal of the
#'   American Statistical Association, 90(431), 1055-1061, Section 6.1.
#' @examples
#' summary(bushfires$area)
"bushfires"

#' Shortest survival times of mice
#'
#' The 84 shortest of the survival times of 208 mice (Kimball 1960), the
#' second example of Wang (1995), who studies the lower tail through the
#' negatives of the times.
#'
#' @format A data frame with 84 rows and one column, `weeks` (survival time
#'   in weeks), in increasing order; the paper prints `-weeks`.
#' @source Wang, J. Z. (1995). Selection of the k largest order statistics
#'   for the domain of attraction of the Gumbel distribution. Journal of the
#'   American Statistical Association, 90(431), 1055-1061, Section 6.2.
#' @examples
#' summary(mice_survival$weeks)
"mice_survival"

#' Best times in the men's 100 m finals of the 1988 and 1992 Olympic Games
#'
#' The best time of each of the 12 fastest athletes in the men's 100 m of the
#' 1988 (Seoul) and 1992 (Barcelona) Olympic Games, the example of Hall &
#' Wang (1999), who estimate the lower end-point of the distribution.
#'
#' @format A data frame with 12 rows: `athlete`, `time` (seconds) and `year`
#'   of the Games, in increasing order of time (names as printed).
#' @source Hall, P., & Wang, J. Z. (1999). Estimating the end-point of a
#'   probability distribution using minimum-distance methods. Bernoulli,
#'   5(1), 177-189, Table 3.
#' @examples
#' hall_wang1999(olympic_100m$time, k = 11, alpha = 0.1, side = "lower",
#'               nsim = 2e4)
"olympic_100m"

#' Flood exceedances of the Wheaton River
#'
#' The 72 exceedances of the threshold 27.50 m^3/s by the flood peaks of the
#' Wheaton River near Carcross, Yukon, 1958-1984, the example of Choulakian &
#' Stephens (2001), printed rounded to one decimal.
#'
#' @format A data frame with 72 rows and one column, `excess` (m^3/s), in the
#'   order printed.
#' @source Choulakian, V., & Stephens, M. A. (2001). Goodness-of-fit tests
#'   for the generalized Pareto distribution. Technometrics, 43(4), 478-484,
#'   Section 3.
#' @examples
#' choulakian_stephens2001(wheaton_river$excess)
"wheaton_river"
