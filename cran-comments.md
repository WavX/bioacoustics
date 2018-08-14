## Test environments
* OS X install (on travis-ci), R release
* ubuntu 14.04 (on travis-ci), R devel and release
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* Fix an OpenMP failure with clang observed on CRAN's machine and properly link the OpenMP library for both C and C++ compilers.
* Add test to detect Solaris OS if CMake version < 3.6
