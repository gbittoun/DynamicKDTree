# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8.11)
project (KDTree)

set(KNOWN_COMPILER FALSE)

if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -std=c++11 -Wall -Werror")
    set(KNOWN_COMPILER TRUE)
elseif ( MSVC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4 /WX")
    set(KNOWN_COMPILER TRUE)
else()
    message(WARNING "You are using an unknown compiler, some features might not be available")
endif()

add_executable (${PROJECT_NAME}_Test main.cpp KDTree.h)

if ( KNOWN_COMPILER )
    add_executable (${PROJECT_NAME}_Test_Optimized main.cpp KDTree.h)

    if ( CMAKE_COMPILER_IS_GNUCC )
        message("Optimizing for GNUCC")
        target_compile_options(${PROJECT_NAME}_Test_Optimized PUBLIC "-O6")
    elseif ( MSVC )
        message("Optimizing for MSVC")
        target_compile_options(${PROJECT_NAME}_Test_Optimized PUBLIC "/O2")
    endif()
else()
    message(WARNING "No known compiler, optimization target will not be generated")
endif()
