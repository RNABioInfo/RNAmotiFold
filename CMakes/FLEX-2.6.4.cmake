include(ExternalProject)

set (FLEX_PREFIX ${CMAKE_BINARY_DIR}/FLEX-prefix)
set (FLEX_LOCAL True)

ExternalProject_Add(
    FLEX
    PREFIX ${FLEX_PREFIX}
    BUILD_IN_SOURCE 1
    DOWNLOAD_EXTRACT_TIMESTAMP true
    URL "https://github.com/westes/flex/releases/download/v2.6.4/flex-2.6.4.tar.gz"
    CONFIGURE_COMMAND ./configure --prefix=${FLEX_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)