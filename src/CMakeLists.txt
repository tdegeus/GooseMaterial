
# define project information
project(GooseSolid)
set(GOOSESOLID_VERSION_NUMBER 0.0.2)

# required to specify the c++ standard
cmake_minimum_required(VERSION 3.0)

# required for install
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# set c++ standard -> c++14
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD LESS 14)
  set(CMAKE_CXX_STANDARD 14)
endif()

# set paths where 'GooseSolid' will be installed
# - root
set(GOOSESOLID_ROOT_DIR ${CMAKE_INSTALL_PREFIX})
# - headers
set(GOOSESOLID_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/${INCLUDE_INSTALL_DIR}")
# - headers
set(INCLUDE_INSTALL_DIR
  "${CMAKE_INSTALL_INCLUDEDIR}/GooseSolid"
  CACHE PATH "Path, relative to CMAKE_PREFIX_PATH, where the GooseSolid header files are installed"
)
# - CMake configuration
set(CMAKEPACKAGE_INSTALL_DIR
  "${CMAKE_INSTALL_DATADIR}/GooseSolid"
  CACHE PATH "Path where GooseSolidConfig.cmake is installed"
)
# - pkg-config
set(PKGCONFIG_INSTALL_DIR
  "${CMAKE_INSTALL_DATADIR}/pkgconfig"
  CACHE PATH "Path where GooseSolid.pc is installed"
)

# list all header files
set(GOOSESOLID_HEADERS
  GooseSolid/ElasticPlasticPotential.h
  GooseSolid/LinearElastic.h
  GooseSolid/LinearElastic_ViscousFluid.h
  GooseSolid/Miscellaneous.h
  GooseSolid/NonLinearElastic.h
  GooseSolid/PlasticLinearElastic.h
  GooseSolid/ViscoPlasticHardeningLinearElastic.h
  GooseSolid/ViscoPlasticLinearElastic.h
)

# Disable pkg-config for native Windows builds
if(NOT WIN32 OR NOT CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
  option(GOOSESOLID_BUILD_PKGCONFIG "Build pkg-config .pc file for GooseSolid" ON)
endif()

# build pkg-config file -> fill in relevant values
if(GOOSESOLID_BUILD_PKGCONFIG)
  configure_file(GooseSolid.pc.in GooseSolid.pc @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GooseSolid.pc DESTINATION ${PKGCONFIG_INSTALL_DIR})
endif()

# configure CMake file
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/GooseSolidConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/GooseSolidConfig.cmake
  PATH_VARS GOOSESOLID_INCLUDE_DIR GOOSESOLID_ROOT_DIR
  INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# install
# - CMake-file
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/GooseSolidConfig.cmake
  DESTINATION
  ${CMAKEPACKAGE_INSTALL_DIR}
)
# - headers
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/${GOOSESOLID_HEADERS}
  DESTINATION
  ${INCLUDE_INSTALL_DIR}/GooseSolid
)

# print information to screen
message(STATUS "")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "|")
message(STATUS "| Use 'make install' to install in '${CMAKE_INSTALL_PREFIX}'")
message(STATUS "| To specify a custom directory call")
message(STATUS "|   cmake . -DCMAKE_INSTALL_PREFIX=yourprefix")
message(STATUS "| For custom paths, add the following line to your '~/.bashrc'")
message(STATUS "|   export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/share/pkgconfig:$PKG_CONFIG_PATH")
message(STATUS "|")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "")
