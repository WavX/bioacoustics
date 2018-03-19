if(getRversion() < "3.3.0")
{
  stop("Your version of R is too old. This package requires R-3.3.0 or newer on Windows.")
}

if (!file.exists("../windows/include/fftw3.h"))
{
  download.file("https://github.com/wavx/rpkg-libs/raw/master/fftw3/fftw3.zip", "lib.zip", quiet = TRUE)

  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows", files = c("include/fftw3.h", paste0("lib/", .Platform$r_arch, "/libfftw3.a")))
  unlink("lib.zip")
}

if (!file.exists("../windows/include/samplerate.h"))
{
  download.file("https://github.com/wavx/rpkg-libs/raw/master/libsamplerate/libsamplerate.zip", "lib.zip", quiet = TRUE)

  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows", files = c("include/samplerate.h", paste0("lib/", .Platform$r_arch, "/libsamplerate.a")))
  unlink("lib.zip")
}

