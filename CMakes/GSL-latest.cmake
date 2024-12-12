include(ExternalProject)

set(GSL_PREFIX ${CMAKE_BINARY_DIR}/GSL-prefix)
set(GSL_LOCAL True)

ExternalProject_Add(
    GSL
    PREFIX ${GSL_PREFIX}
    BUILD_IN_SOURCE 1
    DOWNLOAD_EXTRACT_TIMESTAMP true
    URL "https://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz"
    CONFIGURE_COMMAND ./configure --prefix=${GSL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)