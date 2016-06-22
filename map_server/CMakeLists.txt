cmake_minimum_required(VERSION 2.8.3)
project(map_server)

# GCC specifics
if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -std=c++0x")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb")
endif()

find_package(PkgConfig)
find_package(OpenCV REQUIRED)

include_directories(
  include
  include/map_server
  ${OpenCV_INCLUDE_DIR}
#  ${EIGEN3_INCLUDE_DIR}
)

add_library(map_server
  src/map_server.cpp
)

# test executable
add_executable(test_map_server
  src/map_server.cpp
  examples/test_main.cpp
)
target_link_libraries(test_map_server
  ${OpenCV_LIBS}
)

install(TARGETS map_server
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/obj
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)

install(DIRECTORY include/map_server/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)

install(DIRECTORY examples/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/examples
)
