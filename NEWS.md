# version 0.1.6

* Internal fixes
  * Fix compiler flags
  * Fix magnitude calculation
  * Fix tutorial vignette
  * Update ZC metadata reader for ZC generated with new version of Kaleidoscope

# version 0.1.5

* Internal fixes
  * Fix an OpenMP failure with clang observed on CRAN's machine
  * Properly link the OpenMP library for both C and C++ compilers.
  * Add test to detect Solaris OS if CMake version < 3.6
  
# version 0.1.4

* Link explicitly librt on Solaris
* Improve configure script with better error messages
* Upgrade libsoxr-lsr to 0.1.3.9000

# version 0.1.3

* Add a tutorial (vignette)
* Make the use of OpenMP conditional
* Fix for Solaris

# version 0.1.2

* Automate fftw installation if a suitable version was not found
* Fix installation on OSX
* Fix overloads on Solaris

# version 0.1.1

* Minor bug fixes release

# version 0.1.0

* First release
