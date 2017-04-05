cmake_minimum_required (VERSION 2.6)

enable_language(Fortran)
project (LoopModeling)

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

add_library (loop_closure loop_closure.f90)
target_link_libraries(loop_closure globals geometry mathfunction sturm)

add_library (loop_modeling loop_modeling.f90)
target_link_libraries(loop_modeling globals random mathfunction geometry loop_closure)

add_executable (LoopModeling main_loop_modeling.f90)
target_link_libraries (LoopModeling globals logger random in_out geometry loop_modeling)
