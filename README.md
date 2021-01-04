[![CRAN](http://www.r-pkg.org/badges/version/bioacoustics)](https://cran.r-project.org/package=bioacoustics)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/bioacoustics)](https://cran.r-project.org/package=bioacoustics)
[![R build status](https://github.com/WavX/bioacoustics/workflows/R-CMD-check/badge.svg)](https://github.com/WavX/bioacoustics/actions?workflow=R-CMD-check)

# bioacoustics: detect and extract automatically acoustic features in audio recordings

## Description

bioacoustics contains all the necessary functions to read audio recordings of various formats,
filter noisy files, display audio signals, detect and extract automatically acoustic features
for further analysis such as species identification based on classification of animal vocalizations.

It can be subdivided into three main components:

* Read, extract data (not yet implemented), display, and write Zero-Crossing files.
* Stand-alone tools to convert MP3, WAV, and WAC files.
* Read, display, MP3, WAV or WAC files, filter, and extract automatically acoustic features.


## Installing

* Install stable version from CRAN:
```r
install.packages("bioacoustics")
```

* Install development version from GitHub:
```r
# install.packages("devtools")
devtools::install_github("wavx/bioacoustics", build_vignettes = TRUE)
```

### Windows

Installing bioacoustics from source works under windows when [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is installed. This downloads the system requirements from [rpkg-libs](https://github.com/wavx/rpkg-libs). 

### Linux

For Unix-alikes, [FFTW](http://www.fftw.org/) (>= 3.3.1) is required.

### Contributing

* Contributions are more than welcome, issues and pull requests are the preferred ways of sharing them.

### Authors and contributors

**Authors:** Jean Marchal, Francois Fabianek, Christopher Scott

**Contributors:** Chris Corben, David Riggs, Peter Wilson

### Licence

* **GPLv3**

