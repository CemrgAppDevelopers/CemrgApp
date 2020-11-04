/*=========================================================================

Program:   Medical Imaging & Interaction Toolkit
Module:    $RCSfile$
Language:  C++
Date:      $Date$
Version:   $Revision: 13820 $

Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/ for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include <mitkBaseApplication.h>

#include <QStringList>
#include <QVariant>

int main(int argc, char** argv) {

    mitk::BaseApplication myApp(argc, argv);
    myApp.setSingleMode(true);
    myApp.setApplicationName("CemrgApp v2.1");
    myApp.setOrganizationName("KCL");

    // -------------------------------------------------------------------
    // Here you can switch to your customizable application:
    // -------------------------------------------------------------------

    QStringList preloadLibs;
    preloadLibs << "liborg_mitk_gui_qt_common";
    myApp.setPreloadLibraries(preloadLibs);
    //myApp.setProperty(mitk::BaseApplication::PROP_APPLICATION, "org.mitk.qt.extapplication");
    myApp.setProperty(mitk::BaseApplication::PROP_APPLICATION, "kcl.cemrgapp.mainapp");
    return myApp.run();
}
