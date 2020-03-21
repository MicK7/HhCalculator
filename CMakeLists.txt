cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

set(PluginName HhCalculator)
project(${PluginName})

include(GNUInstallDirs)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(BUILD_SHARED_LIBS ON)

if (NOT ParaView_BINARY_DIR)
  find_package(ParaView REQUIRED)
endif()

# Scan the plugin file in order to set up internal data structures for building
# plugins.
paraview_plugin_scan(
  # The `paraview.plugin` file describing the plugin.
  PLUGIN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Plugin/paraview.plugin"
  # A result variable for the (enabled) plugins found during the scan.
  PROVIDES_PLUGINS plugins
  # Enable plugins during this scan by default.
  ENABLE_BY_DEFAULT ON)

paraview_plugin_build(
  PLUGINS ${plugins})
