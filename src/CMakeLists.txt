cmake_minimum_required (VERSION 2.6)

enable_language(Fortran)
project (LoopModeling)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
    "Choose the installation directory"
    FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
# gfortran
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
# ifort
    set (CMAKE_Fortran_FLAGS_RELEASE "-fast -fpp -heap-arrays -static -no-prec-div -nowarn")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-C -O0 -fpp -heap-arrays -static -no-prec-div -warn all")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
# g77
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

add_library (globals globals.f90)

add_library (logger logger.f90)
target_link_libraries(logger globals)

add_library (random random.f90)
target_link_libraries(random globals)

add_library (mathfunction mathfunction.f90)
target_link_libraries(mathfunction globals)

add_library (geometry geometry.f90)
target_link_libraries(geometry globals logger mathfunction)

add_library (in_out in_out.f90)
target_link_libraries(in_out globals logger geometry)

add_library (sturm sturm.c)

add_library (cluster cluster.f90)
target_link_libraries(cluster globals)

add_library (loop_closure loop_closure.f90)
target_link_libraries(loop_closure globals geometry mathfunction sturm)

add_library (loop_modeling loop_modeling.f90)
target_link_libraries(loop_modeling globals mathfunction geometry loop_closure)

add_executable (LoopModeling main_loop_modeling.f90)
target_link_libraries (LoopModeling globals logger random in_out geometry cluster loop_modeling)

set (EXECUTABLES "LoopModeling")
install (TARGETS ${EXECUTABLES} RUNTIME DESTINATION "bin")
