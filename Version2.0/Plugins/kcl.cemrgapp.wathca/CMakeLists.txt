project(kcl_cemrgapp_wathca)

mitk_create_plugin(
  EXPORT_DIRECTIVE WATHCA_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkSegmentation MitkQtWidgetsExt MitkCemrgAppModule
  PACKAGE_DEPENDS ITK
)
