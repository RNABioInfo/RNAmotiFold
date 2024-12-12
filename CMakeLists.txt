cmake_minimum_required(VERSION 3.16)
project(RNAmotiFold)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#define package lookup and installation macro
macro(test_pack PACKAGE VERSION)
    set(INCLUDE_CMAKE ${CMAKE_CURRENT_SOURCE_DIR}/CMakes/${PACKAGE}-${VERSION}.cmake)
    include(Find${PACKAGE} ${VERSION})
    if(NOT ${${PACKAGE}_FOUND})
        include(${INCLUDE_CMAKE})
    else()
        set (${PACKAGE}_LOCAL False)
    endif()
endmacro()


include(FindGSL)
if (NOT ${GSL_FOUND})
    include(${CMAKE_CURRENT_SOURCE_DIR}/CMakes/GSL-latest.cmake)
endif()

test_pack(BISON 3.8)
test_pack(FLEX 2.6.4)
test_pack(Boost 1.8.6)

include(${CMAKE_CURRENT_SOURCE_DIR}/CMakes/gapcM-latest.cmake)
