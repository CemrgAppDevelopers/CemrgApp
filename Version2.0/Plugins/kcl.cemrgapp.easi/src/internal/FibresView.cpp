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
#include <mitkUnstructuredGrid.h>
#include <mitkImage.h>
#include "EASIView.h"
#include "FibresView.h"


// VTK
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkGlyph3D.h>
#include <vtkSphereSource.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkProperty.h>
#include <vtkCellPicker.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkStreamTracer.h>
#include <vtkTubeFilter.h>
#include <vtkPointSource.h>
#include <vtkCenterOfMass.h>
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
  connect(m_Controls.btnX_checkVtk, SIGNAL(clicked()), this, SLOT(CheckForVtk()));
  connect(m_Controls.btnY_visualise, SIGNAL(clicked()), this, SLOT(Visualise()));
  connect(m_Controls.btnZ_visualiseFibres, SIGNAL(clicked()), this, SLOT(VisualiseFibres()));

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
  meshActor = vtkSmartPointer<vtkActor>::New();
  renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground(0,0,0);

  vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow =
          vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
  m_Controls.widget_1->SetRenderWindow(renderWindow);
  m_Controls.widget_1->GetRenderWindow()->AddRenderer(renderer);

  if (QFile::exists(vtkPath)){
      QMessageBox::warning(NULL, "Attention", "VTK file found. Reading and displaying.");
      PreSurf();
      Visualiser();
      m_Controls.widget_1->GetRenderWindow()->Render();
  } else{
      QMessageBox::warning(NULL, "Attention", "No VTK file yet. Create it by pressing the buttons.");
      MITK_INFO << "VTK file not found. Create it first.";
      m_Controls.btnY_visualise->setEnabled(false);
      m_Controls.btnZ_visualiseFibres->setEnabled(false);
      m_Controls.combo_selector->setEnabled(false);
  }
}

// SLOTS
void FibresView::Preprocessing(){
    bool redoProcess = true;
    if(IsPreprocessingDone()){
        int redoreply = QMessageBox::question(NULL, "Question",
            "Preprocessing files (_shift, _COG, _Cav) found. \nRe-do this step?",
            QMessageBox::Yes, QMessageBox::No);
        if(redoreply == QMessageBox::No){
            redoProcess = false;
        }
    }

    if(redoProcess){
        MITK_INFO <<"[FibresView] Preprocessing Step";
        imagePath = QFileDialog::getOpenFileName(NULL, "Open MASK.nii file",
        dir.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());

        MITK_INFO << "[FibresView] Shifting points to image coordinates.";

        QMessageBox::warning(NULL, "Attention", "Operations may take a few minutes.");
        this->BusyCursorOn();
        mitk::ProgressBar::GetInstance()->AddStepsToDo(2);

        CemrgCommonUtils::OriginalCoordinates(imagePath, ptsFull, shiftPts);
        mitk::ProgressBar::GetInstance()->Progress();

        MITK_INFO << "[FibresView] Centre of Gravity calculation";
        CemrgCommonUtils::CalculateCentreOfGravity(shiftPts, elemFull, cogPts);
        mitk::ProgressBar::GetInstance()->Progress();
        this->BusyCursorOff();

        dilatedCavPath = QFileDialog::getOpenFileName(NULL, "Open dilatedCav.nii file",
        dir.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());

        MITK_INFO << "[FibresView] Region Mapping";
        this->BusyCursorOn();
        CemrgCommonUtils::RegionMapping(dilatedCavPath, cogPts, elemFull, cavElem);
        this->BusyCursorOff();

        FinishedProcess("Preprocessing");
        preprocessDone = true;
    }

    QMessageBox::warning(NULL, "Attention", "Continuing for VTK creation.");
    this->BusyCursorOn();
    CemrgCommonUtils::CarpToVtk(cavElem, shiftPts, vtkPath);
    this->BusyCursorOff();
    FinishedProcess("VTK creation");
}

void FibresView::ExtractSurfaces(){
    bool redoProcess = true;
    if(IsExtractSurfDone()){
        int redoreply = QMessageBox::question(NULL, "Question",
            "Surface files (.surf.vtx) found. \nRe-do this step?",
            QMessageBox::Yes, QMessageBox::No);

        if(redoreply == QMessageBox::No){
            redoProcess = false;
        }

    }

    if(redoProcess){
        ExtractSurfacesProcess();
    }
}

void FibresView::LaplaceSolves(){
    bool redoProcess = true;
    if(IsLaplaceSolvesDone()){
        int redoreply = QMessageBox::question(NULL, "Question",
            "Laplace solves .dat files found. \nRe-do this step?",
            QMessageBox::Yes, QMessageBox::No);

        if(redoreply == QMessageBox::No){
            redoProcess = false;
        }
    }

    if(redoProcess){
        // LaplaceSolvesProcess();
        LaplaceSolvesOpenCarpProcess();
    }

    // laplace solves rectifying to be in the [0,1] range
    QMessageBox::warning(NULL, "Attention", "Checking laplace solve files.");
    CemrgCommonUtils::RectifyFileValues(lapApex);
    CemrgCommonUtils::RectifyFileValues(lapEpi);
    CemrgCommonUtils::RectifyFileValues(lapLV);
    CemrgCommonUtils::RectifyFileValues(lapRV);

    if(laplaceSolvesDone){
        int vtkreply = QMessageBox::question(NULL, "Question",
            "Do you want to add the laplace solves scalar fields to the VTK path?",
        QMessageBox::Yes, QMessageBox::No);
        if(vtkreply == QMessageBox::Yes){
            this->BusyCursorOn();
            mitk::ProgressBar::GetInstance()->AddStepsToDo(4);

            std::vector<double> scalarFieldApex = ReadInScalarField(lapApex);
            CemrgCommonUtils::AppendScalarFieldToVtk(vtkPath, "APEX", "POINT", scalarFieldApex);
            mitk::ProgressBar::GetInstance()->Progress();

            std::vector<double> scalarFieldEpi = ReadInScalarField(lapEpi);
            CemrgCommonUtils::AppendScalarFieldToVtk(vtkPath, "EPI", "POINT", scalarFieldEpi, false);
            mitk::ProgressBar::GetInstance()->Progress();

            std::vector<double> scalarFieldLV = ReadInScalarField(lapLV);
            CemrgCommonUtils::AppendScalarFieldToVtk(vtkPath, "LV", "POINT", scalarFieldLV, false);
            mitk::ProgressBar::GetInstance()->Progress();

            std::vector<double> scalarFieldRV = ReadInScalarField(lapRV);
            CemrgCommonUtils::AppendScalarFieldToVtk(vtkPath, "RV", "POINT", scalarFieldRV, false);
            mitk::ProgressBar::GetInstance()->Progress();
            this->BusyCursorOff();
            FinishedProcess("VTK file Laplace Solves");
        } else{
            QMessageBox::warning(NULL, "Attention", "VTK file not found!");
        }
    }
}

void FibresView::GenerateFibres(){
    bool redoProcess=true;
    if(IsFibresGenerationDone()){
        int redoreply = QMessageBox::question(NULL, "Question",
            "Fibres (_mod.lon) file found. \nRe-do this step?",
            QMessageBox::Yes, QMessageBox::No);

        if(redoreply == QMessageBox::No){
            redoProcess = false;
        }
    }

    if(redoProcess){
        GenerateFibresProcess();
    }

    if (fibresDone){
        int vtkreply = QMessageBox::question(NULL, "Question",
        "Do you want to add the fibres vector field to the VTK path?",
        QMessageBox::Yes, QMessageBox::No);
        if(vtkreply == QMessageBox::Yes){
            this->BusyCursorOn();
            mitk::ProgressBar::GetInstance()->AddStepsToDo(2);
            std::vector<double> fibresVectorField = ReadInVectorField(fibresFile);
            mitk::ProgressBar::GetInstance()->Progress();
            CemrgCommonUtils::AppendVectorFieldToVtk(vtkPath, "fibres", "CELL", fibresVectorField, true);
            mitk::ProgressBar::GetInstance()->Progress();
            this->BusyCursorOff();
            FinishedProcess("VTK file Fibres");
        } else{
            QMessageBox::warning(NULL, "Attention", "Fibres were not added!");
        }
    }

}

// ui slots
void FibresView::CheckForVtk(){
    if(QFile::exists(vtkPath)){
        PreSurf();
        Visualise();

        m_Controls.btnY_visualise->setEnabled(true);
        m_Controls.btnZ_visualiseFibres->setEnabled(true);
        m_Controls.combo_selector->setEnabled(true);
    } else{
        MITK_INFO << "VTK file not found. Create it first.";
        QMessageBox::warning(NULL, "Attention", "Use the buttons above to run the processes and create the vtk file.");
    }
}

void FibresView::Visualise(){
    int comboIndex = m_Controls.combo_selector->currentIndex();
    QString comboText = m_Controls.combo_selector->currentText();
    if(comboIndex==-1){
        comboIndex = 0;
    }
    MITK_INFO << ("[" + QString::number(comboIndex) + "] Combo box text selected: " + comboText ).toStdString();
    ScalarFieldSelector(comboIndex);
    Visualiser(comboIndex);
    m_Controls.widget_1->GetRenderWindow()->Render();

}

void FibresView::VisualiseFibres(){
    if(!fibresVtkCreated){
        SetStreamTracerForFibres();
    }

    MITK_INFO << "Fibres Visualiser";
    vtkSmartPointer<vtkDataSetMapper> meshMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    meshMapper->SetInputData(tubeFilter->GetOutput());

    double max_scalar, min_scalar;
    min_scalar = 0;
    max_scalar = 1;

    meshMapper->SetScalarRange(min_scalar, max_scalar);

    vtkSmartPointer<vtkActor> meshActor = vtkSmartPointer<vtkActor>::New();
    meshActor->SetMapper(meshMapper);
    meshActor->GetProperty()->SetOpacity(1);
    renderer->AddActor(meshActor);

    m_Controls.widget_1->GetRenderWindow()->Render();
    MITK_INFO << "Fibres visualiser finished";
}

// Processes - Actual implementations
void FibresView::ExtractSurfacesProcess(){
    if(IsPreprocessingDone()){
        MITK_INFO << "[FibresView] Creating operations with surface tags.";
        QDialog* inputs = new QDialog(0,0);
        m_UITags.setupUi(inputs);
        connect(m_UITags.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
        connect(m_UITags.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
        int dialogCode = inputs->exec();

        if (dialogCode == QDialog::Accepted) {
            tag_lvendo = m_UITags.txt_lvendo->text();
            tag_lvepi = m_UITags.txt_lvepi->text();
            tag_lvbase = m_UITags.txt_lvbase->text();
            tag_rvendo = m_UITags.txt_rvendo->text();
            tag_rvepi = m_UITags.txt_rvepi->text();
            tag_rvbase = m_UITags.txt_rvbase->text();
            tag_apex = m_UITags.txt_apex->text();

            CheckTag(tag_lvendo, "LVendo");
            CheckTag(tag_lvepi, "LVepi");
            CheckTag(tag_lvbase, "LVbase");
            CheckTag(tag_rvendo, "RVendo");
            CheckTag(tag_rvepi, "RVepi");
            CheckTag(tag_rvbase, "RVbase");
            CheckTag(tag_apex, "Apex");

        } else {
            QMessageBox::warning(NULL, "Attention", "Process cancelled.");
            return;
        }

        // Create operations' strings
        QString opApex, opBIVbase, opLVendo, opRVendo, opEpi;
        QStringList vapex = {tag_apex, tag_lvepi};
        QStringList vbivbase = {tag_lvbase, tag_rvbase, tag_rvepi, tag_lvepi};
        QStringList vlvendo = {tag_lvendo, tag_lvepi, tag_lvbase, tag_apex, tag_rvbase, tag_rvepi, tag_rvendo};
        QStringList vrvendo = {tag_rvendo, tag_lvepi, tag_lvbase, tag_apex, tag_rvbase, tag_rvepi, tag_lvendo};
        QStringList vepi = {tag_lvepi, tag_apex, tag_rvepi, tag_lvbase, tag_rvbase, tag_lvendo, tag_rvendo};

        opApex = CreateSurfExtractOp(vapex, (QStringList() << "-"));
        opBIVbase = CreateSurfExtractOp(vbivbase,  (QStringList() << "," << "-" << ","));
        opLVendo = CreateSurfExtractOp(vlvendo, (QStringList()<<"-"<<","<<","<<","<<","<<","));
        opRVendo = CreateSurfExtractOp(vrvendo, (QStringList()<<"-"<<","<<","<<","<<","<<","));
        opEpi = CreateSurfExtractOp(vepi, (QStringList()<<","<<","<<"-"<<","<<","<<","));

        MITK_INFO << "[FibresView] Extracting BIVbase and Apex";
        // QFile::copy(elemFull, modPathName+".elem");
        // QFile::copy(ptsFull, modPathName+".pts");

        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetDockerImage("alonsojasl/cemrg-meshtool:v1.0");
        MITK_INFO << ("Input: " + name + ", Output:" + modName).toStdString();
        surfApex = cmd->DockerSurfaceFromMesh(dir, name, modName, opApex, "_apex");
        surfBase = cmd->DockerSurfaceFromMesh(dir, name, modName, opBIVbase, "_BIVbase");

        // QFile::remove(modPathName+".elem");
        // QFile::copy(cavElem, modPathName+".elem");
        MITK_INFO << "Copying .pts file to _Cav.pts for other surfaces extraction";
        QString cavName =  name + "_Cav";
        QString cavPathName = dir + mitk::IOUtil::GetDirectorySeparator() + cavName;
        QFile::copy(ptsFull, cavPathName+".pts");

        MITK_INFO << "[FibresView] Extracting LVendo, RVendo and Epi";
        MITK_INFO << ("Input: " + cavName + ", Output:" + modName).toStdString();
        surfLVendo = cmd->DockerSurfaceFromMesh(dir, cavName, modName, opLVendo, "_LVendo");
        surfRVendo = cmd->DockerSurfaceFromMesh(dir, cavName, modName, opRVendo, "_RVendo");
        surfEpi = cmd->DockerSurfaceFromMesh(dir, cavName, modName, opEpi, "_epi");

        FinishedProcess("Extract Surfaces", "Cleaning up extra files");
        extractSurfaceDone = true;

        MITK_INFO << "Cleaning up (.neubc)";
        QFile::remove(modPathName+"_apex.neubc");
        QFile::remove(modPathName+"_BIVbase.neubc");
        QFile::remove(modPathName+"_LVendo.neubc");
        QFile::remove(modPathName+"_RVendo.neubc");
        QFile::remove(modPathName+"_epi.neubc");
        QFile::remove(dir+mitk::IOUtil::GetDirectorySeparator()+name+".fcon");
        QFile::remove(dir+mitk::IOUtil::GetDirectorySeparator()+cavName+".fcon");
        QFile::remove(dir+mitk::IOUtil::GetDirectorySeparator()+cavName+".pts");

        FinishedProcess("Cleaning up .neubc files");

    } else{
        QMessageBox::warning(NULL, "Attention",
            "Preprocessing files (_shift, _COG, _Cav) not found! \nDo the Preprocessing Step");
    }
}

void FibresView::LaplaceSolvesOpenCarpProcess(){
    if(IsExtractSurfDone()){
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetDockerImage("docker.opencarp.org/opencarp/opencarp:latest");

        MITK_INFO << "Copying shift and cav files to _mod";
        QFile::copy(shiftPts, modPathName+".pts");

        if(QFile::exists(cavElem)){
            QFile::remove(modPathName+".elem");
            QFile::rename(cavElem, modPathName+".elem");
        }
        QStringList setToZero, setToOne, regionLabels;
        regionLabels << tag_lvendo << tag_lvepi << tag_lvbase << tag_rvendo << tag_rvepi << tag_rvbase << tag_apex;
        QMessageBox::warning(NULL, "Attention", "Operations may take several minutes.");

        this->BusyCursorOn();
        mitk::ProgressBar::GetInstance()->AddStepsToDo(5);
        mitk::ProgressBar::GetInstance()->Progress();
        MITK_INFO << "[LAPSOLVES] Apex to base";
        setToZero << surfApex.left(surfApex.lastIndexOf(".vtx"));
        setToOne << surfBase.left(surfBase.lastIndexOf(".vtx"));
        QString lapApex_temp = cmd->OpenCarpDockerLaplaceSolves(dir, modName, "lap_apex", setToZero, setToOne, regionLabels);
        setToZero.clear(); setToOne.clear();

        mitk::ProgressBar::GetInstance()->Progress();

        MITK_INFO << "[LAPSOLVES] Epi";
        setToZero << surfLVendo.left(surfLVendo.lastIndexOf(".vtx"));
        setToZero << surfRVendo.left(surfRVendo.lastIndexOf(".vtx"));
        setToOne << surfEpi.left(surfEpi.lastIndexOf(".vtx"));
        QString lapEpi_temp = cmd->OpenCarpDockerLaplaceSolves(dir, modName, "lap_epi", setToZero, setToOne, regionLabels);
        setToZero.clear(); setToOne.clear();

        mitk::ProgressBar::GetInstance()->Progress();

        MITK_INFO << "[LAPSOLVES] LV";
        setToZero << surfEpi.left(surfEpi.lastIndexOf(".vtx"));
        setToZero << surfRVendo.left(surfRVendo.lastIndexOf(".vtx"));
        setToOne << surfLVendo.left(surfLVendo.lastIndexOf(".vtx"));
        QString lapLV_temp = cmd->OpenCarpDockerLaplaceSolves(dir, modName, "lap_lv", setToZero, setToOne, regionLabels);
        setToZero.clear(); setToOne.clear();

        mitk::ProgressBar::GetInstance()->Progress();

        MITK_INFO << "[LAPSOLVES] RV";
        setToZero << surfEpi.left(surfEpi.lastIndexOf(".vtx"));
        setToZero << surfLVendo.left(surfLVendo.lastIndexOf(".vtx"));
        setToOne << surfRVendo.left(surfRVendo.lastIndexOf(".vtx"));
        QString lapRV_temp = cmd->OpenCarpDockerLaplaceSolves(dir, modName, "lap_rv", setToZero, setToOne, regionLabels);
        setToZero.clear(); setToOne.clear();

        mitk::ProgressBar::GetInstance()->Progress();
        this->BusyCursorOff();

        FinishedProcess("Laplace Solves", "Checking for inconsistencies");

        laplaceSolvesDone = true;

    } else{
        QMessageBox::warning(NULL, "Attention",
            "Surface files (.surf.vtx) not found! \nDo the Extract Surface Step");
    }
}

void FibresView::LaplaceSolvesProcess(){
    if(IsExtractSurfDone()){
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetDockerImage("docker.opencarp.org/opencarp/opencarp:latest");

        MITK_INFO << "Copying shift and cav files to _mod";
        QFile::copy(shiftPts, modPathName+".pts");

        if(QFile::exists(cavElem)){
            QFile::remove(modPathName+".elem");
            QFile::rename(cavElem, modPathName+".elem");
        }

        QString paramfile = CemrgCommonUtils::M3dlibLapSolvesParamFile(dir, "param.par", modName, dir, false);

        QMessageBox::warning(NULL, "Attention", "Operations may take several minutes.");
        QString outLapSolves = cmd->ExecuteLaplaceSolves(dir, modName, modName, paramfile);
        FinishedProcess("Laplace Solves", "Cleaning up extra files");

        MITK_INFO << ("[FibresView] Output Laplace Solves: " + outLapSolves).toStdString();

        MITK_INFO(cmd->IsOutputSuccessful(lapApex)) << "[FibresView][LAPSOLVES] Successful Apex";
        MITK_INFO(cmd->IsOutputSuccessful(lapEpi)) << "[FibresView][LAPSOLVES] Successful Epi";
        MITK_INFO(cmd->IsOutputSuccessful(lapLV)) << "[FibresView][LAPSOLVES] Successful LV";
        MITK_INFO(cmd->IsOutputSuccessful(lapRV)) << "[FibresView][LAPSOLVES] Successful RV";

        laplaceSolvesDone = true;
        MITK_INFO << "Cleaning up (.grad)";
        QFile::remove(modPathName+"_lap_apex.grad");
        QFile::remove(modPathName+"_lap_epi.grad");
        QFile::remove(modPathName+"_lap_lv.grad");
        QFile::remove(modPathName+"_lap_rv.grad");
        FinishedProcess("Cleaning Laplace Solves gradient files (.grad)");

    } else{
        QMessageBox::warning(NULL, "Attention",
            "Surface files (.surf.vtx) not found! \nDo the Extract Surface Step");
    }
}

void FibresView::GenerateFibresProcess(){
    if(IsLaplaceSolvesDone()){
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetDockerImage("alonsojasl/cemrg-fibres:v1.0");
        // Default values. Add gui to change them
        MITK_INFO << "[FibresView] Selecting angles values from UI.";
        QDialog* inputs = new QDialog(0,0);
        m_UIAngles.setupUi(inputs);
        connect(m_UIAngles.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
        connect(m_UIAngles.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
        int dialogCode = inputs->exec();

        if (dialogCode == QDialog::Accepted) {
            bool ok1, ok2, ok3, ok4;
            alpha_endo = m_UIAngles.txt_aendo->text().toDouble(&ok1);
            alpha_epi = m_UIAngles.txt_aepi->text().toDouble(&ok2);
            beta_endo = m_UIAngles.txt_bendo->text().toDouble(&ok3);
            beta_epi = m_UIAngles.txt_bepi->text().toDouble(&ok4);

            if(!ok1) SetDefaultAngleAlphaEndo();
            if(!ok2) SetDefaultAngleAlphaEpi();
            if(!ok3) SetDefaultAngleBetaEndo();
            if(!ok4) SetDefaultAngleBetaEpi();
        } else {
            QMessageBox::warning(NULL, "Attention", "Process cancelled.");
            return;
        }

        fibresFile = cmd->DockerComputeFibres(dir, modName, lapApex, lapEpi, lapLV, lapRV, "biv", alpha_endo, alpha_epi, beta_endo, beta_epi);
        FinishedProcess("Fibres generation", "Beginning normalisation of vectors.");

        CemrgCommonUtils::NormaliseFibreFiles(fibresFile, modPathName+"_fibres_norm.lon");
        fibresDone = true;
        FinishedProcess("Fibres normalisation");
    } else{
        QMessageBox::warning(NULL, "Attention", "Laplace Solves files (.dat) not found or have the wrong name.");
    }
}

// helper functions
void FibresView::SetNames(){
    dir = FibresView::directory;
    name = FibresView::basename;
    MITK_INFO << "[FibresView] SetNames()";

    elemFull = dir + mitk::IOUtil::GetDirectorySeparator() + name + ".elem";
    ptsFull =  dir + mitk::IOUtil::GetDirectorySeparator() + name + ".pts";

    std::ifstream fElem(elemFull.toStdString());
    std::ifstream fPts(ptsFull.toStdString());
    fElem >> nElem;
    fPts >> nPts;
    fElem.close();
    fPts.close();

    shiftPts = dir + mitk::IOUtil::GetDirectorySeparator() + name + "_shift.pts";
    cogPts = dir + mitk::IOUtil::GetDirectorySeparator() + name + "_COG.pts";
    cavElem = dir + mitk::IOUtil::GetDirectorySeparator() + name + "_Cav.elem";

    // modPathName will be used in many scenarios
    modName = name + "_mod";
    modPathName = dir + mitk::IOUtil::GetDirectorySeparator() + modName;

    surfApex = modPathName + "_apex.surf.vtx";
    surfBase =  modPathName + "_BIVbase.surf.vtx";
    surfLVendo = modPathName + "_LVendo.surf.vtx";
    surfRVendo = modPathName + "_RVendo.surf.vtx";
    surfEpi = modPathName + "_epi.surf.vtx";

    lapApex = modPathName + "_lap_apex_potential.dat";
    lapEpi = modPathName + "_lap_epi_potential.dat";
    lapLV = modPathName + "_lap_lv_potential.dat";
    lapRV = modPathName + "_lap_rv_potential.dat";

    fibresFile = modPathName + "_fibres.lon";

    vtkPath = modPathName + ".vtk";

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

QString FibresView::CreateSurfExtractOp(QStringList vect, QStringList joinings){
    QString out = "";
    int njoins = joinings.size();
    int nvect = vect.size();
    if (njoins+1 == nvect){
        for(int ix=0; ix<njoins; ix++){
            out += (vect.at(ix) + joinings.at(ix));
        }
        out += vect.at(njoins);
    }
    return out;
}

void FibresView::SetDefaultTag(QString tagname){
    if(tagname.compare("LVendo", Qt::CaseInsensitive) == 0){
        SetDefaultTagLVendo();
    } else if(tagname.compare("LVepi", Qt::CaseInsensitive) == 0){
        SetDefaultTagLVepi();
    } else if(tagname.compare("LVbase", Qt::CaseInsensitive) == 0){
        SetDefaultTagLVbase();
    } else if(tagname.compare("RVendo", Qt::CaseInsensitive) == 0){
        SetDefaultTagRVendo();
    } else if(tagname.compare("RVepi", Qt::CaseInsensitive) == 0){
        SetDefaultTagRVepi();
    } else if(tagname.compare("RVbase", Qt::CaseInsensitive) == 0){
        SetDefaultTagRVbase();
    } if(tagname.compare("Apex", Qt::CaseInsensitive) == 0){
        SetDefaultTagApex();
    }
}

void FibresView::CheckTag(QString someTag, QString tagname){
    bool isOK;
    int numTest = someTag.toInt(&isOK);
    if(!isOK){
        if(!someTag.contains(",") && !someTag.contains("-")){
            MITK_INFO << ("Reverting " + tagname + ": " + someTag + " to default value").toStdString();
            SetDefaultTag(tagname);
        }
    }
}

std::vector<double> FibresView::ReadInScalarField(QString fieldPath){
    std::ifstream fi(fieldPath.toStdString());
    std::vector<double> field(nPts, 0.0);

    for (int ix = 0; ix < nPts; ix++) {
        if(fi.eof()){
            MITK_INFO << "File finished prematurely.";
            break;
        }
        fi >> field[ix];
    }

    return field;
}

std::vector<double> FibresView::ReadInVectorField(QString fieldPath){
    std::ifstream fi(fieldPath.toStdString());
    std::vector<double> field(3*nElem, 0.0);

    int numVect;
    fi >> numVect;

    double discardX, discardY, discardZ;
    for (int ix = 0; ix < nElem; ix++) {
        if(fi.eof()){
            MITK_INFO << "File finished prematurely.";
            break;
        }

        fi >> field[ix + 0*nElem];
        fi >> field[ix + 1*nElem];
        fi >> field[ix + 2*nElem];

        if(numVect > 1){
            fi >> discardX;
            fi >> discardY;
            fi >> discardZ;
        }
    }

    return field;
}

void FibresView::FinishedProcess(QString finishedProcessName, QString xtramessage){
    QString msg = "Finished " + finishedProcessName + " process.";
    msg += (!xtramessage.isEmpty()) ? "" : ("\n" + xtramessage);

    QMessageBox::warning(NULL, "Attention", msg);
}


bool FibresView::IsPreprocessingDone(){
    bool isPreproc;
    isPreproc = QFile::exists(shiftPts) && QFile::exists(cogPts) && QFile::exists(cavElem);
    preprocessDone = isPreproc;

    return isPreproc;
}

bool FibresView::IsExtractSurfDone(){
    bool isExtractSurf;
    bool isSurfFileCorrect = false;

    isExtractSurf = QFile::exists(surfApex) &&
                    QFile::exists(surfBase) &&
                    QFile::exists(surfEpi) &&
                    QFile::exists(surfLVendo) &&
                    QFile::exists(surfRVendo);

    if(isExtractSurf){
        MITK_INFO << "Check: surface files number of points.";
        isSurfFileCorrect = (CemrgCommonUtils::GetTotalFromCarpFile(surfApex)>0) &&
                        (CemrgCommonUtils::GetTotalFromCarpFile(surfBase)>0) &&
                        (CemrgCommonUtils::GetTotalFromCarpFile(surfEpi)>0) &&
                        (CemrgCommonUtils::GetTotalFromCarpFile(surfLVendo)>0) &&
                        (CemrgCommonUtils::GetTotalFromCarpFile(surfRVendo)>0);
    }

    MITK_INFO(!isExtractSurf) << "Surface files not created.";
    MITK_INFO(!isSurfFileCorrect) << "Surface files empty.";

    extractSurfaceDone = (isExtractSurf) ? isSurfFileCorrect : false;

    return extractSurfaceDone;
}

bool FibresView::IsLaplaceSolvesDone(){
    bool isLapSolves;
    isLapSolves = QFile::exists(lapApex) && QFile::exists(lapEpi) && QFile::exists(lapLV) && QFile::exists(lapRV);
    MITK_INFO(!QFile::exists(lapApex)) << "[LAPSOLVES] Unsuccessful Apex";
    MITK_INFO(!QFile::exists(lapEpi)) << "[LAPSOLVES] Unsuccessful Epi";
    MITK_INFO(!QFile::exists(lapLV)) << "[LAPSOLVES] Unsuccessful LV";
    MITK_INFO(!QFile::exists(lapRV)) << "[LAPSOLVES] Unsuccessful RV";

    laplaceSolvesDone = isLapSolves;

    return isLapSolves;
}

bool FibresView::IsFibresGenerationDone(){
    bool isFibres;
    isFibres = QFile::exists(fibresFile);
    fibresDone = isFibres;

    return isFibres;
}

// visualiser functions
void FibresView::ScalarFieldSelector(int vtkIndex){
    MITK_INFO << ("ScalarFieldSelector - index: " + QString::number(vtkIndex)).toStdString();
    int vtkGridSuccess;
    vtkGridSuccess = reader->GetOutput()->GetPointData()->SetActiveScalars(reader->GetScalarsNameInFile(vtkIndex));
    MITK_INFO << reader->GetScalarsNameInFile(vtkIndex);
    MITK_INFO << ("ScalarFieldSelector success index: [" + QString::number(vtkGridSuccess) + "]").toStdString();

}

void FibresView::SetStreamTracerForFibres(){
    MITK_INFO << "Creating stream tracer to display fibres";
    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(5);
    mitk::ProgressBar::GetInstance()->Progress();
    double centre[3];
    vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter =
    vtkSmartPointer<vtkCenterOfMass>::New();
    centerOfMassFilter->SetInputData(reader->GetOutput());
    centerOfMassFilter->SetUseScalarsAsWeights(false);
    centerOfMassFilter->Update();
    centerOfMassFilter->GetCenter(centre);
    mitk::ProgressBar::GetInstance()->Progress();

    std::cout << "COG = (";
    for (int ix = 0; ix < 3; ix++) {
        std::cout << centre[ix] << ",";
    }
    std::cout << ")" << '\n';

    vtkSmartPointer<vtkPointSource> pointSource = vtkSmartPointer<vtkPointSource>::New();
    pointSource->SetCenter(centre);
    pointSource->SetRadius(85000);
    pointSource->SetNumberOfPoints(10000);
    pointSource->Update();
    mitk::ProgressBar::GetInstance()->Progress();

    const char fibreAttribute[] = "fibres";
    reader->GetOutput()->GetPointData()->SetActiveVectors(fibreAttribute);

    vtkSmartPointer<vtkStreamTracer> streamTracer = vtkSmartPointer<vtkStreamTracer>::New();
    streamTracer->SetInterpolatorTypeToDataSetPointLocator();
    streamTracer->SetInputDataObject(reader->GetOutput());
    streamTracer->SetSourceConnection(pointSource->GetOutputPort());
    streamTracer->SetIntegratorTypeToRungeKutta45();
    streamTracer->SetIntegrationDirectionToBoth();
    streamTracer->SetStartPosition(centre);
    streamTracer->Update();
    mitk::ProgressBar::GetInstance()->Progress();

    tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
    tubeFilter->SetRadius(500);
    tubeFilter->SetVaryRadiusToVaryRadiusOff();
    tubeFilter->SetNumberOfSides(3);
    tubeFilter->SetInputConnection(streamTracer->GetOutputPort());
    tubeFilter->Update();

    fibresVtkCreated = true;
    mitk::ProgressBar::GetInstance()->Progress();
    this->BusyCursorOff();
}

// private functions
void FibresView::PreSurf(){
    MITK_INFO << "Reading in VTK file.";

    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(3);
    reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(vtkPath.toStdString().c_str());
    reader->ReadAllScalarsOn();
    reader->ReadAllVectorsOn();
    mitk::ProgressBar::GetInstance()->Progress();
    reader->Update();
    mitk::ProgressBar::GetInstance()->Progress();
    // vtkGrid = reader->GetOutput();

    int numOfScalars = reader->GetNumberOfScalarsInFile();
    int numOfVectors = reader->GetNumberOfVectorsInFile();

    MITK_INFO << ("number of Scalars in read VTK: " + QString::number(numOfScalars)).toStdString();
    for (int ix = 0; ix < numOfScalars; ix++) {
        // Set combo combo box
        m_Controls.combo_selector->addItem(reader->GetScalarsNameInFile(ix));
        MITK_INFO << reader->GetScalarsNameInFile(ix);
    }

    MITK_INFO << ("number of Vectors in read VTK: " + QString::number(numOfVectors)).toStdString();
    for (int ix = 0; ix < numOfVectors; ix++) {
        MITK_INFO << reader->GetVectorsNameInFile(ix);
    }

    fibresVtkCreated = false;
    mitk::ProgressBar::GetInstance()->Progress();
    this->BusyCursorOff();
}

void FibresView::Visualiser(int vtkIndex){
    MITK_INFO << ("Visualiser - index: " + QString::number(vtkIndex));
    vtkSmartPointer<vtkDataSetMapper> meshMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    meshMapper->SetInputData(reader->GetOutput());

    double max_scalar, min_scalar;
    if(vtkIndex==0){
        meshMapper->SetScalarModeToUseCellData();
        max_scalar = 30;
        min_scalar = 1;
    } else{
        meshMapper->SetScalarModeToUsePointData();
        // meshMapper->ScalarVisibilityOn();
        min_scalar = 0;
        max_scalar = 1;
    }
    meshMapper->SetScalarRange(min_scalar, max_scalar);

    vtkSmartPointer<vtkActor> meshActor = vtkSmartPointer<vtkActor>::New();
    meshActor->SetMapper(meshMapper);
    meshActor->GetProperty()->SetOpacity(1);
    renderer->AddActor(meshActor);
}
