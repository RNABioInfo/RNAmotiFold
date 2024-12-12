include(ExternalProject)
set(Boost_USE_MULTITHREADED ON)
set(Boost_DEBUG OFF)
set(Boost_ROOT ${CMAKE_BINARY_DIR}/boost-prefix)

set (BOOST_LOCAL True)

ExternalProject_Add(
    Boost
    PREFIX ${Boost_ROOT}
    BUILD_IN_SOURCE 1
    DOWNLOAD_EXTRACT_TIMESTAMP true
    URL "https://github.com/boostorg/boost/releases/download/boost-1.86.0/boost-1.86.0-cmake.tar.gz"
    CONFIGURE_COMMAND ./bootstrap.sh
    --prefix=<INSTALL_DIR>
    --with-libraries=test
    --with-libraries=program_options
    --with-libraries=system
    --with-libraries=filesystem
    --with-libraries=serialization
    BUILD_COMMAND ./b2 install link=static variant=release threading=multi runtime-link=static
    INSTALL_COMMAND ""
    INSTALL_DIR ${Boost_ROOT}

)