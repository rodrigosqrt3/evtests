# evtests 0.1.0

First release.

* Estimators of the extreme value index: `hill1975()`, `pickands1975()`,
  `dekkers1989()`, `drees2004()`; generalized Pareto maximum likelihood by
  `grimshaw1993()`.
* Tests of the extreme value condition: `dietrich2002()` (`E` and `PE`, with
  the critical values of Husler & Li 2006) and `drees2006()`, with limiting
  distributions simulated and stored in the package.
* Tests of the Gumbel domain: `hasofer_wang1992()`, `segers_teugels2000()`,
  `neves2006()`, `neves_fraga_alves2007()`, `fraga_alves2016()`,
  `castillo2014()` and the CV-plot `castillo2014_cv()`.
* Finite end-point: `neves_pereira2010()`, `hall_wang1999()`, `smith1985()`
  and `fraga_alves_neves2014()` with the upper bound of Fraga Alves, Neves &
  Rosario (2016).
* Goodness of fit of the generalized Pareto distribution:
  `choulakian_stephens2001()`.
* Choice of the number of excesses or of the threshold: `wang1995()`,
  `drees_kaufmann1998()`, `danielsson2001()`, `guillou_hall2001()`,
  `bader2018()`.
* Data sets from the source papers: `us_cities_1940`, `nylon_fibres`,
  `bushfires`, `mice_survival`, `olympic_100m`, `wheaton_river`.
* Vignette `reproduction`: what was reproduced from each paper and where
  printed values differ from the papers' formulas or designs.
