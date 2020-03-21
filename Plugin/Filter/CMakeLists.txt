set(classes
  vtkHhCalculator
)

set(SERVER_MANAGER_XML_files
  HhCalculator.xml
)

set(SERVER_MANAGER_SOURCES_files
  vtkHhCalculator.cxx
)

set(SERVER_MANAGER_HEADERS_files
  vtkHhCalculator.h
)

vtk_module_add_module(HhCalculator::vtkHhCalculator
  SOURCES ${SERVER_MANAGER_SOURCES_files}
  HEADERS ${SERVER_MANAGER_HEADERS_files}
)

paraview_add_server_manager_xmls(
  XMLS ${SERVER_MANAGER_XML_files}
)

