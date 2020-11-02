/*=========================================================================

Program:   Medical Imaging & Interaction Toolkit
Language:  C++
Date:      $Date$
Version:   $Revision$

Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
/*=========================================================================
 *
 * CemrgApp Common Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbench.h>
#include <berryIWorkbenchWindow.h>
#include <berryIWorkbenchPage.h>

// Qmitk
#include <mitkCoreObjectFactory.h>
#include <mitkProgressBar.h>
#include <QmitkIOUtil.h>

// CemrgAppModule
#include <CemrgCommonUtils.h>
#include "QmitkCemrgAppCommonTools.h"

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QSignalMapper>

const std::string QmitkCemrgAppCommonTools::VIEW_ID = "org.mitk.views.cemrgappcommontools";

void QmitkCemrgAppCommonTools::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::LoadMesh);
    connect(m_Controls.button_2, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::ConvertToCarto);
    connect(m_Controls.button_3, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::ConvertCarpToVtk);
}

void QmitkCemrgAppCommonTools::SetFocus() {
}

void QmitkCemrgAppCommonTools::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void QmitkCemrgAppCommonTools::LoadMesh() {

    QString path = "";
    path = QFileDialog::getOpenFileName(
                NULL, "Open Mesh Data File", QmitkIOUtil::GetFileOpenFilterString());
    CemrgCommonUtils::AddToStorage(
                CemrgCommonUtils::LoadVTKMesh(path.toStdString()), "Mesh", this->GetDataStorage());
}

void QmitkCemrgAppCommonTools::ConvertToCarto() {

    QString path = "";
    path = QFileDialog::getOpenFileName(
                NULL, "Open Mesh Data File", QmitkIOUtil::GetFileOpenFilterString());

    if (path.isEmpty() || !path.endsWith(".vtk")) {
        QMessageBox::warning(NULL, "Attention", "Select Correct Input File!");
        return;
    }

    CemrgCommonUtils::ConvertToCarto(path.toStdString());
    QMessageBox::information(NULL, "Attention", "Conversion Completed!");
}

void QmitkCemrgAppCommonTools::ConvertCarpToVtk(){
    QString pathElem = "";
    QString pathPts = "";
    pathElem = QFileDialog::getOpenFileName(NULL, "Open Mesh .elem File");
    if (pathElem.isEmpty() || !pathElem.endsWith(".elem")) {
        QMessageBox::warning(NULL, "Attention", "Select Correct Input (.elem) File!");
        return;
    }
    QFileInfo fi(pathElem);
    QString dir = fi.absolutePath();
    QString vtkPath = dir + mitk::IOUtil::GetDirectorySeparator() +  fi.baseName() + ".vtk";

    pathPts = QFileDialog::getOpenFileName(NULL, "Open Mesh .pts File", dir.toStdString().c_str());

    if (pathPts.isEmpty() || !pathPts.endsWith(".pts")) {
        QMessageBox::warning(NULL, "Attention", "Select Correct Input (.pts) File!");
        return;
    }

    int regionScalarsReply = QMessageBox::question(NULL, "Question",
            "Include region as (cell) scalar field?", QMessageBox::Yes, QMessageBox::No);

    CemrgCommonUtils::CarpToVtk(pathElem, pathPts, vtkPath, (regionScalarsReply==QMessageBox::Yes));

    int appendScalarFieldReply = QMessageBox::question(NULL, "Question",
            "Append a scalar field from a file?", QMessageBox::Yes, QMessageBox::No);

    if (appendScalarFieldReply==QMessageBox::Yes){
        QString path="";
        QString typeData="";
        int nElem = CemrgCommonUtils::GetTotalFromCarpFile(pathElem);
        int nPts = CemrgCommonUtils::GetTotalFromCarpFile(pathPts);
        int nField;

        while (appendScalarFieldReply==QMessageBox::Yes){
            path = QFileDialog::getOpenFileName(NULL, "Open Scalar field (.dat) file", dir.toStdString().c_str());
            QFileInfo fi2(path);
            std::vector<double> field = CemrgCommonUtils::ReadScalarField(path);

            nField = field.size();
            MITK_INFO << ("FieldSize: " + QString::number(nField)).toStdString();
            if(nField==nElem){
                typeData = "CELL";
            } else if(nField==nPts){
                typeData = "POINT";
            } else {
                MITK_INFO << "Inconsistent file size";
                break;
            }
            CemrgCommonUtils::AppendScalarFieldToVtk(vtkPath, fi2.baseName(), typeData, field);

            appendScalarFieldReply = QMessageBox::question(NULL, "Question",
                    "Append another scalar field from a file?", QMessageBox::Yes, QMessageBox::No);
        }
    }
}
