## Test environments
* Local 
	* R version 4.0.3
	* Platform: x86_64-apple-darwin17.0 (64-bit)
	* Running under: macOS High Sierra 10.13.4
* Travis-CI
	* R versions 3.6.3, 4.0.2, and Under development (unstable) (2021-05-19 r80339)
	* Platform: x86_64-pc-linux-gnu
	* Running under: Ubuntu 16.04.6 LTS
* Appveyor
	* R version 4.0.2 (2020-06-22)
	* Platform: x86_64-w64-mingw32/x64 (64-bit)
	* Running under: Windows Server 2012 R2 x64 (build 9600)
* RHub 
	* Fedora Linux, R-devel, clang, gfortran
	* Ubuntu Linux 20.04 LTS, R-release, GCC
	* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Win-builder
	* R versions 4.0.5, 4.1.0, 4.2.0
	* Platform: x86_64-w64-mingw32 (64-bit) 

RHub check on Windows Server 2008 R2 SP1, R-devel, 32/64 bit failed. 
	* appears to be due to native issues in RHub (Bioconductor yet not
	supporting R-devel).

All checks reporting OK. 

## Downstream dependencies
Nothing to report.

## Additional Notes
Submission of updated package, version 1.1.1. 