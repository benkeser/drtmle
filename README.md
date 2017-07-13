# R/`drtmle`

[![Travis-CI Build Status](https://travis-ci.org/benkeser/drtmle.svg?branch=master)](https://travis-ci.org/benkeser/drtmle)
[![AppVeyor Build  Status](https://ci.appveyor.com/api/projects/status/github/benkeser/drtmle?branch=master&svg=true)](https://ci.appveyor.com/project/benkeser/drtmle)
[![Coverage Status](https://img.shields.io/codecov/c/github/benkeser/drtmle/master.svg)](https://codecov.io/github/benkeser/drtmle?branch=master)
[![CRAN](http://www.r-pkg.org/badges/version/drtmle)](http://www.r-pkg.org/pkg/drtmle)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Doubly Robust, Asymptotically Linear Targeted Minimum Loss-Based Estimation 

---

## Description

`drtmle` is an R package that computes marginal effect estimators for binary 
treatments on continuous and binary outcomes. The TMLE estimators are doubly
robust, not only with respect to consistency, but also with respect to asymptotic 
normality, under assumptions included in [Benkeser et al (2016, submitted)](http://biostats.bepress.com/ucbbiostat/paper356/). 

---

## Installation

- Install the most recent _stable release_:
  `devtools::install_github("benkeser/drtmle")`

- To contribute, install the _development version_:
  `devtools::install_github("benkeser/drtmle", ref = "develop")`

---

## License

&copy; 2016-2017 [David C. Benkeser](http://www.benkeserstatistics.com)

The contents of this repository are distributed under the MIT license. See
below for details:
```
The MIT License (MIT)

Copyright (c) 2016-2017 David C. Benkeser

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```