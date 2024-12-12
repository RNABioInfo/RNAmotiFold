include(ExternalProject)

set (GAPC_PREFIX ${CMAKE_BINARY_DIR}/gapc-prefix)

if (${FLEX_EXTERNAL})
    set (FLEX_PATH ${FLEX_PREFIX}/bin/flex)
else()
    set (FLEX_PATH ${FLEX_EXECUTEABLE})
endif()

if (${BISON_EXTERNAL})
    set (BISON_PATH ${BISON_PREFIX}/bin/bison)
else()
    set (BISON_PATH ${BISON_EXECUTEABLE})
endif()

if (${GSL_EXTERNAL})
    set (GSL_CONFIG ${GSL_PREFIX}/bin/gsl-config)
    set (GSL_PATH ${GSL_PREFIX}/lib/libgsl.so)
else()
    set (GSL_CONFIG ${GSL_CONFIG_EXECUTEABLE})
    set (GSL_PATH ${GSL_LIBRARY})
endif()

if (${BOOST_EXTERNAL})
    set (BOOST_PATH ${Boost_ROOT})
else()
    set (BOOST_PATH ${BOOST_INCLUDE_DIRS})
endif()


ExternalProject_Add(
    gapcM
    PREFIX ${GAPC_PREFIX}
    BUILD_IN_SOURCE 1
    DOWNLOAD_EXTRACT_TIMESTAMP true
    GIT_REPOSITORY "https://github.com/RNABioInfo/gapcM.git"
    CONFIGURE_COMMAND ./configure --prefix=${GAPC_PREFIX} FLEX=${FLEX_PATH} BISON=${BISON_PATH} GSL_CONFIG=${GSL_CONFIG} GSL=${GSL_PATH} --with-boost=${BOOST_PATH}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)
add_dependencies(gapcM GSL BISON FLEX Boost)