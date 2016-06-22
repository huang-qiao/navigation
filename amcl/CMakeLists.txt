cmake_minimum_required(VERSION 2.8.3)
project(amcl)

# GCC specifics
if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -std=c++0x")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb")
endif()

find_package(Eigen3 REQUIRED)

include_directories(
  include
  include/amcl
  include/amcl/map
  include/amcl/sensors
  include/amcl/pf
  ${EIGEN3_INCLUDE_DIR}
)

add_library(amcl_pf
  src/amcl/pf/pf.c
  src/amcl/pf/pf_kdtree.c
  src/amcl/pf/pf_pdf.c
  src/amcl/pf/pf_vector.c
  src/amcl/pf/eig3.c
  src/amcl/pf/pf_draw.c
)

add_library(amcl_map
  src/amcl/map/map.c
  src/amcl/map/map_cspace.cpp
  src/amcl/map/map_range.c
  src/amcl/map/map_store.c
  src/amcl/map/map_draw.c
)

add_library(amcl_sensors
  src/amcl/sensors/amcl_sensor.cpp
  src/amcl/sensors/amcl_odom.cpp
  src/amcl/sensors/amcl_laser.cpp
)

add_library(amcl
  src/amcl/amcl_imp.cpp
)

target_link_libraries(amcl_sensors amcl_map amcl_pf)

target_link_libraries(amcl
  amcl_sensors
  amcl_map
  amcl_pf
)

install(TARGETS amcl amcl_sensors amcl_map amcl_pf
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/obj
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)

install(DIRECTORY include/amcl/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)

install(DIRECTORY examples/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/examples
)
