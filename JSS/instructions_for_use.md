# Reproducing Doubly-Robust Inference in R using `drtmle`

See the rendered version of these instructions at
https://github.com/benkeser/drtmle/blob/jss/instructions_for_use.md.

To reproduce the submitted manuscript in its entirety, we assume the user has 
access to a machine that includes: 

- `R`, GNU Make, git;
- an installation of the `R` package `renv`.

The `drtmle` package and relevant manuscript code can be downloaded from GitHub
as follows.

```bash
# clone the repository
git clone https://github.com/benkeser/drtmle/
# move into the repository
cd drtmle
# commit associated with submission
git checkout 4331441d059a317c2cea61cc1714f5962428af8e
```

Next, use `renv` to synchronize the `R` environment for rendering the
manuscript.

```bash
make jss_environment
```

Finally, the manuscript can be rendered as follows.

```bash
make jss
# ...when finished open the rendered pdf
open jss/using_drtmle_jss.pdf
```

__Note__: If GNU make is not available, the above recipes can be executed 
at the command line as follows.

```bash
# from the repository root
# restore environment
cd JSS & Rscript -e "renv::restore()"
# build manuscript
Rscript -e "rmarkdown::render('using_drtmle_jss.Rmd')"
```