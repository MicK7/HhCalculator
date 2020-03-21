# create a paraview plugin containing server manager xml and the server
paraview_add_plugin(${PluginName}
  REQUIRED_ON_SERVER
  VERSION "1.1"
  MODULES HhCalculator::vtkHhCalculator
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Filter/vtk.module")
