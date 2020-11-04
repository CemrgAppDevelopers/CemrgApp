set(SRC_CPP_FILES
)

set(INTERNAL_CPP_FILES
  kcl_cemrgapp_mainapp_Activator.cpp
  QmitkCemrgApplication.cpp
  QmitkCemrgAppCommonTools.cpp
  QmitkCemrgWorkbenchAdvisor.cpp
  perspectives/QmitkCemrgJBPerspective.cpp
  perspectives/QmitkCemrgHCPerspective.cpp
  perspectives/QmitkCemrgRRPerspective.cpp
  perspectives/QmitkCemrgEasiPerspective.cpp
  perspectives/QmitkCemrgPowertransPerspective.cpp
  perspectives/QmitkCemrgWathcaPerspective.cpp
)

set(UI_FILES
  src/internal/QmitkCemrgAppCartoExport.ui
  src/internal/QmitkCemrgAppCommonToolsControls.ui
)

set(MOC_H_FILES
  src/internal/kcl_cemrgapp_mainapp_Activator.h
  src/internal/QmitkCemrgApplication.h
  src/internal/QmitkCemrgAppCommonTools.h
  src/internal/QmitkCemrgWorkbenchAdvisor.h
  src/internal/perspectives/QmitkCemrgJBPerspective.h
  src/internal/perspectives/QmitkCemrgHCPerspective.h
  src/internal/perspectives/QmitkCemrgRRPerspective.h
  src/internal/perspectives/QmitkCemrgEasiPerspective.h
  src/internal/perspectives/QmitkCemrgPowertransPerspective.h
  src/internal/perspectives/QmitkCemrgWathcaPerspective.h
)

# list of resource files which can be used by the plug-in
# system without loading the plug-ins shared library,
# for example the icon used in the menu and tabs for the
# plug-in views in the workbench
set(CACHED_RESOURCE_FILES
  plugin.xml
  resources/icon.xpm
  resources/icon_research.xpm
)

# list of Qt .qrc files which contain additional resources
# specific to this plugin
set(QRC_FILES
  resources/CemrgApp.qrc
)

set(CPP_FILES )

foreach(file ${SRC_CPP_FILES})
  set(CPP_FILES ${CPP_FILES} src/${file})
endforeach(file ${SRC_CPP_FILES})

foreach(file ${INTERNAL_CPP_FILES})
  set(CPP_FILES ${CPP_FILES} src/internal/${file})
endforeach(file ${INTERNAL_CPP_FILES})
