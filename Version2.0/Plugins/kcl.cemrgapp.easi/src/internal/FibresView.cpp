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

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchPage.h>
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

// Qmitk
#include <QmitkIOUtil.h>
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkNodePredicateProperty.h>
#include <mitkImage.h>
#include "EASIView.h"
#include "FibresView.h"


// VTK
#include <vtkGlyph3D.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkProperty.h>
#include <vtkCellPicker.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCell.h>
#include <vtkMath.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkDecimatePro.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkCellDataToPointData.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>

// Qt
#include <QMessageBox>
#include <QDesktopWidget>
#include <QDir>
#include <QDirIterator>
#include <QFileInfo>
#include <QFileDialog>
#include <QStringList>

// MitkCemrgAppModule
#include <CemrgCommandLine.h>
#include <CemrgCommonUtils.h>

QString FibresView::basename;
QString FibresView::directory;

const std::string FibresView::VIEW_ID = "org.mitk.views.fibres";

// FibresView::~FibresView(){
// }

void FibresView::SetDirectoryFile(const QString directory, const QString basename){
    FibresView::basename = basename;
    FibresView::directory = directory;
}

void FibresView::SetFocus() {
  m_Controls.btn1_preproc->setFocus();
}

void FibresView::OnSelectionChanged(
    berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void FibresView::CreateQtPartControl(QWidget *parent) {
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi(parent);
  connect(m_Controls.btn1_preproc, SIGNAL(clicked()), this, SLOT(Preprocessing()));
  connect(m_Controls.btn2_extractsurf, SIGNAL(clicked()), this, SLOT(ExtractSurfaces()));
  connect(m_Controls.btn3_lapsolves, SIGNAL(clicked()), this, SLOT(LaplaceSolves()));
  connect(m_Controls.btn4_fibres, SIGNAL(clicked()), this, SLOT(GenerateFibres()));
  connect(m_Controls.btnX_cancel, SIGNAL(clicked()), this, SLOT(CancelFibres()));
  connect(m_Controls.btnY_visualise, SIGNAL(clicked()), this, SLOT(Visualise()));

  // Setup names and paths
  SetNames();

  // Set default tags
  SetDefaultTagLVendo();
  SetDefaultTagLVepi();
  SetDefaultTagLVbase();
  SetDefaultTagRVendo();
  SetDefaultTagRVepi();
  SetDefaultTagRVbase();
  SetDefaultTagApex();

  // Set default angles
  SetDefaultAngleAlphaEndo();
  SetDefaultAngleAlphaEpi();
  SetDefaultAngleBetaEndo();
  SetDefaultAngleBetaEpi();

  //Setup renderer
  surfActor = vtkSmartPointer<vtkActor>::New();
  renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground(0,0,0);

  vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow =
          vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
  m_Controls.widget_1->SetRenderWindow(renderWindow);
  m_Controls.widget_1->GetRenderWindow()->AddRenderer(renderer);

  //Setup keyboard interactor
  // callBack = vtkSmartPointer<vtkCallbackCommand>::New();
  // callBack->SetCallback(KeyCallBackFunc);
  // callBack->SetClientData(this);
  // interactor = m_Controls.widget_1->GetRenderWindow()->GetInteractor();
  // interactor->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());
  // interactor->GetInteractorStyle()->KeyPressActivationOff();
  // interactor->GetInteractorStyle()->AddObserver(vtkCommand::KeyPressEvent, callBack);
}

// functionality slots
void FibresView::Preprocessing(){
    SetNames();
    MITK_INFO <<"[FibresView] Preprocessing Step";
    imagePath = QFileDialog::getOpenFileName(NULL, "Open MASK.nii file",
                    dir.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());

    MITK_INFO << "[FibresView] Shifting points to image coordinates.";

    QMessageBox::warning(NULL, "Attention", "Operations may take a few minutes.");
    CemrgCommonUtils::OriginalCoordinates(imagePath, ptsFull, shiftPts);

    MITK_INFO << "[FibresView] Centre of Gravity calculation";
    CemrgCommonUtils::CalculateCentreOfGravity(shiftPts, elemFull, cogPts);

    dilatedCavPath = QFileDialog::getOpenFileName(NULL, "Open dilatedCav.nii file",
                    dir.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());

    MITK_INFO << "[FibresView] Region Mapping";
    CemrgCommonUtils::RegionMapping(dilatedCavPath, cogPts, elemFull, cavElem);

    QMessageBox::warning(NULL, "Attention", "Finished preprocessing Step!");
    preprocessDone = true;
}

void FibresView::ExtractSurfaces(){
    if(!preprocessDone){
        int reply = QMessageBox::question(NULL, "Question",
                    "Do you have the files for preprocessing?",
                    QMessageBox::Yes, QMessageBox::No);
        if(reply == QMessageBox::Yes){
            preprocessDone = true;
            ExtractSurfaces();
        } else{
            QMessageBox::warning(NULL, "Attention", "Do the Preprocessing step!");
        }
    } else {
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetDockerImage("alonsojasl/cemrg-meshtool:v1.0");

        MITK_INFO << "[FibresView] Creating operations with surface tags.";
        QDialog* inputs = new QDialog(0,0);
        m_UITags.setupUi(inputs);
        connect(m_UITags.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
        connect(m_UITags.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
        int dialogCode = inputs->exec();

        if (dialogCode == QDialog::Accepted) {
            bool ok1, ok2, ok3, ok4, ok5, ok6, ok7;
            tag_lvendo = m_UITags.txt_lvendo->text().toInt(&ok1);
            if(!ok1) SetDefaultTagLVendo();
            tag_lvepi = m_UITags.txt_lvepi->text().toInt(&ok2);
            if(!ok2) SetDefaultTagLVepi();
            tag_lvbase = m_UITags.txt_lvbase->text().toInt(&ok3);
            if(!ok3) SetDefaultTagLVbase();
            tag_rvendo = m_UITags.txt_rvendo->text().toInt(&ok4);
            if(!ok4) SetDefaultTagRVendo();
            tag_rvepi = m_UITags.txt_rvepi->text().toInt(&ok5);
            if(!ok5) SetDefaultTagRVepi();
            tag_rvbase = m_UITags.txt_rvbase->text().toInt(&ok6);
            if(!ok6) SetDefaultTagRVbase();
            tag_apex = m_UITags.txt_apex->text().toInt(&ok7);
            if(!ok7) SetDefaultTagApex();
        }

        // Create operations' strings
        QString opApex, opBIVbase, opLVendo, opRVendo, opEpi;
        std::vector<int> vapex = {tag_apex, tag_lvendo};
        std::vector<int> vbivbase = {tag_lvbase, t tag_rvepi, tag_lvepi};
        std::vector<int> vlvendo = {tag_lvendo, tag_lvepi, tag_lvbase, tag_apex, tag_rvbase, tag_rvepi, tag_rvendo};
        std::vector<int> vrvendo = {tag_rvendo, tag_lvepi, tag_lvbase, tag_apex, tag_rvbase, tag_rvepi, tag_lvendo};
        std::vector<int> vepi = {tag_lvepi, tag_apex, tag_rvepi, tag_lvbase, tag_rvbase, tag_lvendo, tag_rvendo};

        opApex = CreateSurfExtractOp(vapex, (QStringList() << "-"));
        opBIVbase = CreateSurfExtractOp(vbivbase,  (QStringList() << "," << "-" << ","));
        opLVendo = CreateSurfExtractOp(vlvendo, (QStringList()<<"-"<<","<<","<<","<<","<<","));
        opRVendo = CreateSurfExtractOp(vrvendo, (QStringList()<<"-"<<","<<","<<","<<","<<","));
        opEpi = CreateSurfExtractOp(vepi, (QStringList()<<","<<","<<"-"<<","<<","<<","));

        MITK_INFO << opApex.toStdString();
        MITK_INFO << opBIVbase.toStdString();
        MITK_INFO << opLVendo.toStdString();
        MITK_INFO << opRVendo.toStdString();
        MITK_INFO << opEpi.toStdString();

        MITK_INFO << "[FibresView] Extracting BIVbase and Apex";
        QFile::copy(elemFull, modPathName+".elem");
        QFile::copy(ptsFull, modPathName+".pts");

        surfApex = cmd->DockerSurfaceFromMesh(dir, modName, "3-1", "_apex");
        surfBase = cmd->DockerSurfaceFromMesh(dir, modName, "2,4-5,1", "_BIVbase");

        QFile::remove(modPathName+".elem");
        QFile::copy(cavElem, modPathName+".elem");

        MITK_INFO << "[FibresView] Extracting LVendo, RVendo and Epi";
        surfLVendo = cmd->DockerSurfaceFromMesh(dir, modName, "10-1,2,3,4,5,30", "_LVendo");
        surfRVendo = cmd->DockerSurfaceFromMesh(dir, modName, "30-1,2,3,4,5,10", "_RVendo");
        surfEpi = cmd->DockerSurfaceFromMesh(dir, modName, "1,3,5-2,4,10,30", "_epi");

        extractSurfaceDone = true;
    }
}

void FibresView::LaplaceSolves(){
    if(!extractSurfaceDone){
        int reply = QMessageBox::question(NULL, "Question",
                    "Do you have the .dat files from the Laplace solves?",
                    QMessageBox::Yes, QMessageBox::No);
        if(reply == QMessageBox::Yes){
            extractSurfaceDone = true;
            LaplaceSolves();
        } else{
            QMessageBox::warning(NULL, "Attention", "Do the Laplace Solves step!");
        }
    } else{
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetDockerImage("alonsojasl/meshtools3d-lapsolves:v1.0");

        MITK_INFO << "Copying shift and cav files to _Mod";
        QFile::copy(shiftPts, modPathName+".pts");

        if(QFile::exists(cavElem)){
            QFile::remove(modPathName+".elem");
            QFile::rename(cavElem, modPathName+".elem");
        }

        QString paramfile = CemrgCommonUtils::M3dlibLapSolvesParamFile(dir, "param.par", modName, dir, false);

        QMessageBox::warning(NULL, "Attention", "Operations may take several minutes.");
        QString outLapSolves = cmd->ExecuteLaplaceSolves(dir, modName, name, paramfile);

        MITK_INFO << ("[FibresView] Output Laplace Solves: " + outLapSolves).toStdString();

        MITK_INFO(cmd->IsOutputSuccessful(lapApex)) << "[FibresView][LAPSOLVES] Successful Apex";
        MITK_INFO(cmd->IsOutputSuccessful(lapEpi)) << "[FibresView][LAPSOLVES] Successful Epi";
        MITK_INFO(cmd->IsOutputSuccessful(lapLV)) << "[FibresView][LAPSOLVES] Successful LV";
        MITK_INFO(cmd->IsOutputSuccessful(lapRV)) << "[FibresView][LAPSOLVES] Successful RV";

        laplaceSolvesDone = true;
    }
}

void FibresView::GenerateFibres(){
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetDockerImage("alonsojasl/cemrg-fibres:v1.0");

    bool lapsolveapex = cmd->IsOutputSuccessful(lapApex);
    bool lapsolveepi = cmd->IsOutputSuccessful(lapEpi);
    bool lapsolvelv = cmd->IsOutputSuccessful(lapLV);
    bool lapsolverv = cmd->IsOutputSuccessful(lapRV);

    if(lapsolveapex && lapsolveepi && lapsolvelv && lapsolverv){
        // Default values. Add gui to change them
        MITK_INFO << "[FibresView] Selecting angles values from UI.";
        QDialog* inputs = new QDialog(0,0);
        m_UIAngles.setupUi(inputs);
        connect(m_UITags.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
        connect(m_UITags.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
        int dialogCode = inputs->exec();

        if (dialogCode == QDialog::Accepted) {
            bool ok1, ok2, ok3, ok4;
            alpha_endo = m_UITags.txt_aendo->text().toDouble(&ok1);
            if(!ok1) SetDefaultAngleAlphaEndo();
            alpha_epi = m_UITags.txt_aepi->text().toDouble(&ok2);
            if(!ok2) SetDefaultAngleAlphaEpi();
            beta_endo = m_UITags.txt_bendo->text().toDouble(&ok4);
            if(!ok3) SetDefaultAngleBetaEndo();
            beta_epi = m_UITags.txt_bepi->text().toDouble(&ok3);
            if(!ok4) SetDefaultAngleBetaEpi();
        }

        cmd->DockerComputeFibres(dir, modName, lapApex, lapEpi, lapLV, lapRV, "biv", alpha_endo, alpha_epi, beta_endo, beta_epi);
    } else{
        MITK_INFO(!lapsolveapex) << "[LAPSOLVES] Unsuccessful Apex";
        MITK_INFO(!lapsolveepi) << "[LAPSOLVES] Unsuccessful Epi";
        MITK_INFO(!lapsolvelv) << "[LAPSOLVES] Unsuccessful LV";
        MITK_INFO(!lapsolverv) << "[LAPSOLVES] Unsuccessful RV";
    }

}

// ui slots
void FibresView::CancelFibres(){

}

void FibresView::Visualise(){

}


// helper functions
void FibresView::SetNames(){
    dir = FibresView::directory;
    name = FibresView::basename;
    MITK_INFO << "[FibresView] SetNames()";

    elemFull = dir + mitk::IOUtil::GetDirectorySeparator() + name + ".elem";
    ptsFull =  dir + mitk::IOUtil::GetDirectorySeparator() + name + ".pts";

    shiftPts = dir + mitk::IOUtil::GetDirectorySeparator() + name + "_shift.pts";
    cogPts = dir + mitk::IOUtil::GetDirectorySeparator() + name + "_COG.pts";
    cavElem = dir + mitk::IOUtil::GetDirectorySeparator() + name + "_Cav.elem";

    // modPathName will be used in many scenarios
    modPathName = dir + mitk::IOUtil::GetDirectorySeparator() + name + "_Mod";
    modName = name + "_Mod";

    lapApex = dir + mitk::IOUtil::GetDirectorySeparator() + name + "_lap_apex_potential.dat";
    lapEpi = dir + mitk::IOUtil::GetDirectorySeparator() + name + "_lap_epi_potential.dat";
    lapLV = dir + mitk::IOUtil::GetDirectorySeparator() + name + "_lap_lv_potential.dat";
    lapRV = dir + mitk::IOUtil::GetDirectorySeparator() + name + "_lap_rv_potential.dat";

    QFile::copy(dir + mitk::IOUtil::GetDirectorySeparator() + name + ".lon", modPathName+".lon");

    preprocessDone = false;
    extractSurfaceDone = false;
    laplaceSolvesDone = false;
    fibresDone = false;

    MITK_INFO << ("[...] Elements path: " + elemFull).toStdString();
    MITK_INFO << ("[...] Points path: " + ptsFull).toStdString();
    MITK_INFO << ("[...] Shift path: " + shiftPts).toStdString();
    MITK_INFO << ("[...] COG path: " + cogPts).toStdString();
    MITK_INFO << ("[...] Cav path: " + cavElem).toStdString();
}

QString FibresView::CreateSurfExtractOp(std::vector<int> vect, QStringList joinings){
    QString out = "";
    int njoins = joinings.size();
    int nvect = vect.size();
    if (njoins+1 == nvect){
        for(int ix=0; ix<njoins; ix++){
            out += (QString::number(vect.at(ix)) + joinings.at(ix));
        }
        out += QString::number(vect.at(njoins));
    }
    return out;
}
