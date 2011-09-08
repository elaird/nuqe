cmake_minimum_required(VERSION 2.6)

project(foo) 

#options here
#set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -Wno-long-long -march=native")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-long-long -march=native")

#get ROOT stuff here
exec_program(root-config ARGS "--incdir" OUTPUT_VARIABLE ROOT_INCLUDE_DIR )
include_directories(${ROOT_INCLUDE_DIR})
#exec_program(root-config ARGS "--libdir" OUTPUT_VARIABLE ROOT_LIBRARY_DIR )
#link_directories(${ROOT_LIBRARY_DIR})
exec_program(root-config ARGS "--libs" OUTPUT_VARIABLE ROOT_LIBS )
#message(STATUS ${ROOT_LIBS})

#project here
set(foo_SRCS src/foo.C src/TT_params.C src/TT_generator.C src/TT_nucleus.C src/TT_drawer.C src/TT_event.C)
add_executable(foo ${foo_SRCS})
target_link_libraries(foo ${ROOT_LIBS})
