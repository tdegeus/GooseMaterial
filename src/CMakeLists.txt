
# define project information
project(GooseMaterial)
set(GOOSEMATERIAL_VERSION_NUMBER 0.0.2)

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

# set paths where 'GooseMaterial' will be installed
# - root
set(GOOSEMATERIAL_ROOT_DIR ${CMAKE_INSTALL_PREFIX})
# - headers
set(GOOSEMATERIAL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/${INCLUDE_INSTALL_DIR}")
# - headers
set(INCLUDE_INSTALL_DIR
  "${CMAKE_INSTALL_INCLUDEDIR}/GooseMaterial"
  CACHE PATH "Path, relative to CMAKE_PREFIX_PATH, where the GooseMaterial header files are installed"
)
# - CMake configuration
set(CMAKEPACKAGE_INSTALL_DIR
  "${CMAKE_INSTALL_DATADIR}/GooseMaterial"
  CACHE PATH "Path where GooseMaterialConfig.cmake is installed"
)
# - pkg-config
set(PKGCONFIG_INSTALL_DIR
  "${CMAKE_INSTALL_DATADIR}/pkgconfig"
  CACHE PATH "Path where GooseMaterial.pc is installed"
)

# list all header files
set(GOOSEMATERIAL_HEADERS
  GooseMaterial/ElasticPlasticPotential.h
  GooseMaterial/LinearElastic.h
  GooseMaterial/LinearElastic_ViscousFluid.h
  GooseMaterial/Miscellaneous.h
  GooseMaterial/NonLinearElastic.h
  GooseMaterial/PlasticLinearElastic.h
  GooseMaterial/ViscoPlasticHardeningLinearElastic.h
  GooseMaterial/ViscoPlasticLinearElastic.h
)

# Disable pkg-config for native Windows builds
if(NOT WIN32 OR NOT CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
  option(GOOSEMATERIAL_BUILD_PKGCONFIG "Build pkg-config .pc file for GooseMaterial" ON)
endif()

# build pkg-config file -> fill in relevant values
if(GOOSEMATERIAL_BUILD_PKGCONFIG)
  configure_file(GooseMaterial.pc.in GooseMaterial.pc @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GooseMaterial.pc DESTINATION ${PKGCONFIG_INSTALL_DIR})
endif()

# configure CMake file
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/GooseMaterialConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/GooseMaterialConfig.cmake
  PATH_VARS GOOSEMATERIAL_INCLUDE_DIR GOOSEMATERIAL_ROOT_DIR
  INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# install
# - CMake-file
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/GooseMaterialConfig.cmake
  DESTINATION
  ${CMAKEPACKAGE_INSTALL_DIR}
)
# - headers
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/${GOOSEMATERIAL_HEADERS}
  DESTINATION
  ${INCLUDE_INSTALL_DIR}/GooseMaterial
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
