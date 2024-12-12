cmake_minimum_required(VERSION 3.16)
project(RNAmotiFold)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#define package lookup and installation macro
macro(test_pack PACKAGE REQ_VERSION AVAIL_VERSION)
    set(CMAKE_EXTERNAL_PACKAGE ${CMAKE_CURRENT_SOURCE_DIR}/CMakes/${PACKAGE}-${AVAIL_VERSION}.cmake)
    set(CMAKE_LOCAL_PACKAGE ${CMAKE_CURRENT_SOURCE_DIR}/CMakes/${PACKAGE}_local.cmake)
    find_package(${PACKAGE} ${REQ_VERSION})
    if(NOT ${${PACKAGE}_FOUND})
        include(${CMAKE_EXTERNAL_PACKAGE})
    else()
        set (${PACKAGE}_EXTERNAL False)
        add_executable(${PACKAGE} IMPORTED GLOBAL)
    endif()
endmacro()

macro(test_pack_boost PACKAGE REQ_VERSION AVAIL_VERSION BOOST_LIBS)
    set(CMAKE_EXTERNAL_PACKAGE ${CMAKE_CURRENT_SOURCE_DIR}/CMakes/${PACKAGE}-${AVAIL_VERSION}.cmake)
    set(CMAKE_LOCAL_PACKAGE ${CMAKE_CURRENT_SOURCE_DIR}/CMakes/${PACKAGE}_local.cmake)
    find_package(${PACKAGE} ${REQ_VERSION} COMPONENTS ${BOOST_LIBS})
    if(NOT ${${PACKAGE}_FOUND})
        include(${CMAKE_EXTERNAL_PACKAGE})
    else()
        set (${PACKAGE}_EXTERNAL False)
        add_executable(${PACKAGE} IMPORTED GLOBAL)
    endif()
endmacro()


test_pack(BISON 2.4.1 3.8)
test_pack(FLEX 2.5.34 2.6.4)
test_pack(GSL 2.7 2.8)
test_pack_boost(Boost 1.36 1.8.6 list(test,program_options,system,filesystem,serialization))

include(${CMAKE_CURRENT_SOURCE_DIR}/CMakes/gapcM-latest.cmake)
