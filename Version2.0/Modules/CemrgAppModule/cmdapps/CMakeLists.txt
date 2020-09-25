option(BUILD_CemrgCMDApps "Build ALL command-line apps for Cemrg Module" OFF)
option(BUILD_CEMRG_FIXSHELL "Build fixing shell values command line app. " ON)
option(BUILD_CEMRG_CLIPPINGTOOL "Build clipping valve command line app. " ON)
option(BUILD_CEMRG_PROJECTLGE "Build lge projection command line app. " ON)
option(BUILD_CEMRG_IM2INR "Build image to inr command line app. " ON)
option(BUILD_CEMRG_VENTRICLE_SEGMENTATION_RELABEL "Build ventricle segmentation relabelling command line app" ON)
option(BUILD_CEMRG_CARP_UTILS "Build CemrgApp+Carp utils command line app" ON)
option(BUILD_CEMRG_LABELS_TO_SURF "Build CemrgApp extract surface from labels utils command line app" ON)

if(BUILD_CemrgCMDApps)
  mitkFunctionCreateCommandLineApp(
    NAME CemrgTemplateApp
    DEPENDS MitkCemrgAppModule
    CPP_FILES CemrgAppCmdTemplate.cpp
  )
endif()

if(BUILD_CEMRG_FIXSHELL)
  mitkFunctionCreateCommandLineApp(
    NAME CemrgFixShell
    DEPENDS MitkCemrgAppModule
    CPP_FILES CemrgFixShellApp.cpp
  )
endif()

if(BUILD_CEMRG_CLIPPINGTOOL)
  mitkFunctionCreateCommandLineApp(
    NAME CemrgClippingTool
    DEPENDS MitkCemrgAppModule
    CPP_FILES CemrgClippingTool.cpp
  )
endif()

if(BUILD_CEMRG_PROJECTLGE)
  mitkFunctionCreateCommandLineApp(
    NAME ProjectLgeToVtkMesh
    DEPENDS MitkCemrgAppModule
    CPP_FILES CemrgProjectLgeToVtkMesh.cpp
  )
endif()

if(BUILD_CEMRG_IM2INR)
  mitkFunctionCreateCommandLineApp(
    NAME CemrgIM2INR
    DEPENDS MitkCemrgAppModule
    CPP_FILES CemrgIM2INR.cpp
  )
endif()

if(BUILD_CEMRG_VENTRICLE_SEGMENTATION_RELABEL)
  mitkFunctionCreateCommandLineApp(
    NAME CemrgVentricleSegmRelabel
    DEPENDS MitkCemrgAppModule
    CPP_FILES CemrgVentricleSegmRelabel.cpp
  )
endif()

if(BUILD_CEMRG_CARP_UTILS)
    mitkFunctionCreateCommandLineApp(
      NAME CemrgCarpUtils
      DEPENDS MitkCemrgAppModule
      CPP_FILES CemrgCarpUtils.cpp
    )
endif()

if(BUILD_CEMRG_LABELS_TO_SURF)
    mitkFunctionCreateCommandLineApp(
        NAME CemrgVtx2Surf
        DEPENDS MitkCemrgAppModule
        CPP_FILES CemrgVtxToSurf.cpp
    )
endif()
