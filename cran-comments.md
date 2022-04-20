## Test environments
* Local 
	* R version 4.0.3
	* Platform: x86_64-apple-darwin17.0 (64-bit)
	* Running under: macOS High Sierra 10.13.4 and macOS Big Sur 11.6
* RHub 
	* Fedora Linux, R-devel, clang, gfortran
	* Ubuntu Linux 20.04 LTS, R-release, GCC
	* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Win-builder
	* devel, release, oldrelease
	* R versions 4.0.5, 4.1.0, 4.2.0
	* Platform: x86_64-w64-mingw32 (64-bit) 

Most checks reporting OK.

Some checks report 2 NOTES that appear safe to ignore (copied below):

>  Maintainer: 'David Benkeser <benkeser@emory.edu>'
>  
>  New submission
>  
>  Package was archived on CRAN
>  
>  Possibly misspelled words in DESCRIPTION:
>    Benkeser (13:27)
>    Laan (14:13)
>    MJ (13:79)
>    Nonparametric (2:22)
>    al (13:39)
>    der (14:9)
>    et (13:36)
>  
>  CRAN repository db overrides:
>    X-CRAN-Comment: Archived on 2022-04-04 as 'coercion to logical'
>      errors were not corrected in time.

## Downstream dependencies
Nothing to report.

## Additional Notes
Submission of updated package, version 1.1.1.

The package was archived on 2022-04-04 due to the 'coercion to logical' error.
As far as I can tell from the output, this error actually occurs in the {earth}
package. To avoid the build error, I have removed this package from the vignette
examples and from the "Suggests" field of DESCRIPTION.