# Builds the data sets bundled with evtests, each transcribed from the
# published source cited in R/data.R.
# Run from the package root: source("data-raw/datasets.R")

# Hill (1975), Table 1 (p. 1173): the 30 largest cities of the United States
# in 1940, in thousands of inhabitants, from the Statistical Abstract of the
# United States, 1950, p. 57.
us_cities_1940 <- data.frame(
  rank = 1:30,
  population = 1000 * c(7455, 3397, 1931, 1623, 1504, 878, 859, 816, 771, 672,
                        663, 635, 587, 576, 495, 492, 456, 430, 399, 387, 385,
                        368, 325, 322, 319, 306, 305, 302.3, 302.2, 301.2)
)

# Grimshaw (1993), Table 1 (p. 186): exceedances (kg/mm^2) of a testing
# threshold by the tensile strength of 15 nylon carpet fibres, in the order
# printed (three rows of five).
nylon_fibres <- data.frame(
  excess = c(.051, .140, .365, .561, .030,
             .268, .184, .100, .876, .092,
             .011, .200, .518, .338, .056)
)

# Wang (1995), Section 6.1 (p. 1059): areas (hectares) of the 30 largest of
# 75 Australian bushfires recorded in 1986-1987, from Smith (1993).
bushfires <- data.frame(
  area = c(5800, 3500, 2300, 2160, 1000, 900, 800, 700, 450, 400,
           360, 300, 300, 300, 200, 200, 200, 200, 150, 150,
           120, 110, 100, 99, 82, 80, 80, 55, 50, 40)
)

# Wang (1995), Section 6.2 (p. 1060): the 84 shortest survival times (weeks)
# of 208 mice (Kimball 1960), printed as negatives in decreasing order.
mice_survival <- data.frame(
  weeks = c(40, 48, 50, 54, 56, 59, 62, 63, 67, 67, 69, 70,
            71, 73, 73, 76, 77, 80, 81, 81, 82, 83, 84, 86,
            86, 87, 88, 88, 88, 88, 88, 89, 90, 90, 91, 93,
            94, 95, 96, 97, 97, 98, 99, 99, 100, 100, 100, 100,
            101, 101, 101, 102, 102, 103, 103, 103, 103, 103, 104, 104,
            104, 105, 105, 106, 106, 106, 107, 108, 109, 109, 110, 110,
            110, 111, 111, 111, 112, 113, 113, 114, 114, 115, 116, 116)
)

# Hall & Wang (1999), Table 3 (p. 185): best time of each of the 12 fastest
# athletes in the men's 100 m finals of the 1988 and 1992 Olympic Games.
olympic_100m <- data.frame(
  athlete = c("C. Lewis", "L. Christie", "C. Smith", "F. Fredericks",
              "D. Mitchell", "B. Surin", "L. Burrell", "R. da Silva",
              "D. Williams", "O. Adenniken", "R. Stewart", "D. Ezinwa"),
  time = c(9.92, 9.96, 9.99, 10.02, 10.04, 10.09, 10.10, 10.11, 10.11, 10.12,
           10.22, 10.26),
  year = c(1988, 1992, 1988, 1992, 1992, 1992, 1992, 1988, 1988, 1992, 1992,
           1992)
)

# Choulakian & Stephens (2001), Section 3 (p. 481): the 72 exceedances
# (m^3/s) of the threshold 27.50 by the flood peaks of the Wheaton River near
# Carcross, Yukon, 1958-1984, in the order printed.
wheaton_river <- data.frame(
  excess = c(1.7, 2.2, 14.4, 1.1, 0.4, 20.6, 5.3, 0.7, 1.9, 13, 12, 9.3, 1.4,
             18.7, 8.5, 25.5, 11.6, 14.1, 22.1, 1.1, 2.5, 14.4, 1.7, 37.6,
             0.6, 2.2, 39, 0.3, 15, 11, 7.3, 22.9, 1.7, 0.1, 1.1, 0.6, 9, 1.7,
             7, 20.1, 0.4, 2.8, 14.1, 9.9, 10.4, 10.7, 30, 3.6, 5.6, 30.8,
             13.3, 4.2, 25.5, 3.4, 11.9, 21.5, 27.6, 36.4, 2.7, 64, 1.5, 2.5,
             27.4, 1, 27.1, 20.2, 16.8, 5.3, 9.7, 27.5, 2.5, 27)
)

save(us_cities_1940, file = "data/us_cities_1940.rda", compress = "xz")
save(wheaton_river, file = "data/wheaton_river.rda", compress = "xz")
save(olympic_100m, file = "data/olympic_100m.rda", compress = "xz")
save(nylon_fibres, file = "data/nylon_fibres.rda", compress = "xz")
save(bushfires, file = "data/bushfires.rda", compress = "xz")
save(mice_survival, file = "data/mice_survival.rda", compress = "xz")
