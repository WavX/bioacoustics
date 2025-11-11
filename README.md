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

## System Requirements

This package requires the following system dependencies:
- **C++11**: A C++11 compatible compiler
- **FFTW3**: Fast Fourier Transform library (>= 3.3.1)
- **GNU Make**: Build automation tool

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

### System Dependencies Installation

#### Windows
Installing bioacoustics from source works under windows when [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is installed. This downloads the system requirements from [rpkg-libs](https://github.com/wavx/rpkg-libs).

#### Linux
For Unix-alikes, [FFTW](http://www.fftw.org/) (>= 3.3.1) is required.

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install build-essential libfftw3-dev
```

**CentOS/RHEL:**
```bash
sudo yum install gcc-c++ fftw-devel
```

**Fedora:**
```bash
sudo dnf install gcc-c++ fftw-devel
```

#### macOS
```bash
# Using Homebrew
brew install fftw

# Using MacPorts
sudo port install fftw-3
```

### Contributing
* Contributions are more than welcome, issues and pull requests are the preferred ways of sharing them.

### Authors and contributors
**Authors:** Jean Marchal, Francois Fabianek, Christopher Scott

**Contributors:** Chris Corben, David Riggs, Peter Wilson

### Licence
* **GPLv3**
