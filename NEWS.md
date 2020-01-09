# drtmle 1.0.5

As of January 2019:
* Version 1.0.5 released on GitHub and CRAN.
* Addition of `snow` and `data.table` packages to suggests to quell re-occurring
   warnings and errors on CRAN builds.
* Added `nSL` option for fitting and averaging multiple super learners as part
   of the estimation procedure. 
* Added `adapt_g` option for outcome-adaptive propensity score fitting. 

# drtmle 1.0.4

As of August 2019:
* Added minor touch-ups and link fixes to documentation and vignettes.
* Improved how slots are ordered upon being included in the return object.

As of July 2019:
* Version 1.0.4.9001 released on GitHub.
* Removed dependency on `plyr` package.
* Added option to average over repeated Super Learner fits.

As of December 18, 2018:
* Version 1.0.4 released on GitHub and CRAN.
* Resolved issues arising from `returnModels` option when users input nuisance
   parameters.
* Added option to bypass `future` parallelization calls for easier debugging.
* Fixed bugs in standard TMLE implementation -- namely, more robust fluctuations
   and corrected variance estimators.

# drtmle 1.0.3

As of July 2, 2018:
* Version 1.0.3 released on GitHub and CRAN.
* Fixed warnings on CRAN builds.

# drtmle 1.0.2

As of February 5, 2018:
* Version 1.0.2 released on GitHub and CRAN.
* Replaced `foreach` parallelization with `future`.
* Included more robust Super Learner methods.
* Fixed test to pass build with long doubles removed.
* Accommodated returning estimated influence functions with `drtmle()` fit for
   power users.
* Incorporated minor documentation corrections and updates.

As of December 11, 2017:
* Version 1.0.2.9000 released on GitHub.
* More robust convex combination SuperLearner implemented.

# drtmle 1.0.0

As of August 17, 2017:
* Version 1.0.0 released on CRAN.
* Version 1.0.0.9000 released on GitHub.

As of August 15, 2017:
* Version 1.0.0 ready for CRAN release.

# drtmle 0.0.1

As of April 05, 2017:
* The first public release of this package (v. 0.0.1) is made available on
   GitHub.
