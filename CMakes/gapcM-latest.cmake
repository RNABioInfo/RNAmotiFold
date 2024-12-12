include(ExternalProject)

set (GAPC_PREFIX ${CMAKE_BINARY_DIR}/gapc-prefix)

if (${FLEX_LOCAL})
    set (FLEX ${FLEX_PREFIX}/bin/flex)
else()
    set (FLEX ${FLEX_EXECUTEABLE})
endif()

if (${BISON_LOCAL})
    set (BISON ${BISON_PREFIX}/bin/bison)
else()
    set (BISON ${BISON_EXECUTEABLE})
endif()

if (${GSL_LOCAL})
    set (GSL ${GSL_PREFIX}/bin/gsl-config)
else()
    set (GSL ${GSL_LIBRARIES})
endif()

if (${BOOST_LOCAL})
    set (BOOST ${Boost_ROOT})
else()
    set (BOOST ${BOOST_INCLUDE_DIRS})
endif()


ExternalProject_Add(
    gapcM
    PREFIX ${GAPC_PREFIX}
    BUILD_IN_SOURCE 1
    DOWNLOAD_EXTRACT_TIMESTAMP true
    GIT_REPOSITORY "https://github.com/RNABioInfo/gapcM.git"
    CONFIGURE_COMMAND ./configure --prefix=${GAPC_PREFIX} FLEX=${FLEX} BISON=${BISON} GSL_CONFIG=${GSL} GSL=${GSL} --with-boost=${BOOST}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)
add_dependencies(gapcM GSL BISON FLEX Boost)