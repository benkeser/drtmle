md:
	Rscript -e "rmarkdown::render('README.Rmd', output_file = 'README.md')"

site:
	Rscript -e "rmarkdown::render('README.Rmd', output_file = 'README.md')"
	Rscript -e "pkgdown::build_site()"

check:
	Rscript -e "devtools::check()"

check_local_rhub:
	Rscript -e "rhub::local_check_linux(getwd(), image = 'rhub/debian-clang-devel')"
	Rscript -e "rhub::local_check_linux(getwd(), image = 'rhub/debian-gcc-devel')"
	Rscript -e "rhub::local_check_linux(getwd(), image = 'rhub/debian-gcc-devel-nold')"
	Rscript -e "rhub::local_check_linux(getwd(), image = 'rhub/debian-gcc-patched')"
	Rscript -e "rhub::local_check_linux(getwd(), image = 'rhub/debian-gcc-release')"
	Rscript -e "rhub::local_check_linux(getwd(), image = 'rhub/fedora-clang-devel')"
	Rscript -e "rhub::local_check_linux(getwd(), image = 'rhub/fedora-gcc-devel')"
	Rscript -e "rhub::local_check_linux(getwd(), image = 'rhub/centos6-epel')"
	Rscript -e "rhub::local_check_linux(getwd(), image = 'rhub/centos6-epel-rdt')"
	Rscript -e "rhub::local_check_linux(getwd(), image = 'rhub/rocker-gcc-san')"
	Rscript -e "rhub::local_check_linux(getwd(), image = 'rhub/ubuntu-gcc-devel')"
	Rscript -e "rhub::local_check_linux(getwd(), image = 'rhub/ubuntu-gcc-release')"
	Rscript -e "rhub::local_check_linux(getwd(), image = 'rhub/ubuntu-rchk')"

checkfast:
	Rscript -e "devtools::check(build_args = '--no-build-vignettes')"

test:
	Rscript -e "devtools::test()"

doc:
	Rscript -e "devtools::document()"

build:
	Rscript -e "devtools::build()"

buildfast:
	Rscript -e "devtools::build(vignettes = FALSE)"

style:
	Rscript -e "styler::style_pkg()"

jss_environment: 
	cd JSS && Rscript -e "renv::restore()"

jss:
	cd JSS && Rscript -e "rmarkdown::render('using_drtmle_jss.Rmd')"

.PHONY: jss
