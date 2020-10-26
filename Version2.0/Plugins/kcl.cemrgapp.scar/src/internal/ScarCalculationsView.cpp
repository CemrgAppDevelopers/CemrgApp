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
 * Scar Quantification
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * jose.solislemus@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchPage.h>
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

// Qmitk
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkNodePredicateProperty.h>
#include <mitkImage.h>
#include "ScarCalculationsView.h"
#include "AtrialScarView.h"

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
#include <QStringList>

// MitkCemrgAppModule
#include <CemrgCommandLine.h>

QString ScarCalculationsView::fileName;
QString ScarCalculationsView::directory;
QString ScarCalculationsView::predir;
QString ScarCalculationsView::postdir;
QString ScarCalculationsView::advdir;
QString ScarCalculationsView::preScarFile;
QString ScarCalculationsView::postScarFile;

const std::string ScarCalculationsView::VIEW_ID = "org.mitk.views.scarcalculations";

ScarCalculationsView::~ScarCalculationsView() {
    //inputs->deleteLater();
}

void ScarCalculationsView::SetCalculationsPaths(const QString directory) {

    ScarCalculationsView::directory = directory;
    ScarCalculationsView::predir = directory +
            mitk::IOUtil::GetDirectorySeparator() + "PRE" +
            mitk::IOUtil::GetDirectorySeparator() + "ANALYSIS";
    ScarCalculationsView::postdir = directory +
            mitk::IOUtil::GetDirectorySeparator() + "POST" +
            mitk::IOUtil::GetDirectorySeparator() + "ANALYSIS";
    ScarCalculationsView::advdir = directory +
            mitk::IOUtil::GetDirectorySeparator() + "ADVANCED_ANALYSIS";
    ScarCalculationsView::preScarFile = "";
    ScarCalculationsView::postScarFile = "";
}

bool ScarCalculationsView::CheckForRequiredFiles() {

    QString searchPre = ScarCalculationsView::predir + mitk::IOUtil::GetDirectorySeparator();
    QString searchPost = ScarCalculationsView::postdir + mitk::IOUtil::GetDirectorySeparator();

    int responsePre = ScarCalculationsView::SearchDirectory(searchPre);
    int responsePost = ScarCalculationsView::SearchDirectory(searchPost);

    bool ret = false;
    if (responsePre + responsePost >= 8){
        ret = true;
    }
    return ret;
}

int ScarCalculationsView::SearchDirectory(QString searchDir){
    MITK_INFO << ("[INFO] Searching files on directory: " + searchDir).toStdString();
    int response = 0;
    bool debugVar = false;
    bool isPre = searchDir.contains("PRE", Qt::CaseSensitive);

    QDirIterator qiter(searchDir, QDirIterator::Subdirectories);
    while(qiter.hasNext()) { // look for .nii LGE and MRA files in pre
        QFileInfo finfo(qiter.next());
        if (finfo.fileName().contains(".nii", Qt::CaseSensitive)) {
            if (finfo.fileName().contains("LGE", Qt::CaseSensitive)){
                MITK_INFO(debugVar) << "[DEBUG] found: LGE.";
                response++;
            }

            if (finfo.fileName().contains("MRA", Qt::CaseSensitive)){
                MITK_INFO(debugVar) << "[DEBUG] found: MRA.";
                response++;
            }
        }
        if (finfo.fileName().contains("prodThresholds", Qt::CaseSensitive)){
            MITK_INFO(debugVar) << "[DEBUG] found: Thresholds file.";
            response++;
        }

        if (finfo.fileName().contains(".vtk", Qt::CaseSensitive)){
            if (!finfo.fileName().contains("Normalised", Qt::CaseSensitive)){
                if (finfo.fileName().contains("MaxScar", Qt::CaseSensitive)){
                    MITK_INFO(debugVar) << "[DEBUG] found: Scar VTK.";
                    if(isPre){
                        preScarFile = finfo.fileName();
                    } else{
                        postScarFile = finfo.fileName();
                    }
                    response++;
                }
            }
        }
    }
    return response;
}

QStringList ScarCalculationsView::CheckForAdvancedDirectoryFiles() {

    // Only checks for MaxScarPre/Post and prodThresholdsPre/Post
    QStringList need2load = {"Pre", "Post"};
    std::vector<int> v;
    QString searchDir = ScarCalculationsView::advdir + mitk::IOUtil::GetDirectorySeparator();
    QDirIterator itdir(searchDir, QDirIterator::Subdirectories);
    while(itdir.hasNext()) { // look for .nii LGE and MRA files in pre
        QFileInfo finfo(itdir.next());

        if (finfo.fileName().contains("prodThresholdsPre.txt", Qt::CaseSensitive))
            need2load.removeAt(0);

        if (finfo.fileName().contains("prodThresholdsPost.txt", Qt::CaseSensitive)) {
            int rmat = ((need2load.count()==2) ? 1 : 0);
            need2load.removeAt(rmat);
        }
    }
    return need2load;
}

void ScarCalculationsView::GetInputsFromFile() {

    MITK_INFO << "GET INPUTS FROM FILE.\n";
    double data1[5]; //, data2[5];
    QString prodPath = QString();
    QString prodPathOut = ScarCalculationsView::advdir + mitk::IOUtil::GetDirectorySeparator();
    QDir advd(ScarCalculationsView::advdir);

    if (advd.mkdir(ScarCalculationsView::advdir)) {
        QMessageBox::warning(NULL, "Advanced analysis folder",
                             ("The advanced analysis folder was created in:\n\n"+
                              ScarCalculationsView::advdir.toStdString()).c_str());
        MITK_INFO << "The advanced analysis folder was created.";
    }

    // QStringList alreadyInFolder = {"Pre", "Post"};
    QStringList need2load = CheckForAdvancedDirectoryFiles();

    if (need2load.count() == 0) {
        MITK_INFO << "Using thresholds alread in ADVANCED_ANALYSIS folder." ;
    }
    else {
        MITK_INFO << "Loading..." + need2load.join(", ").toStdString();

        for (int i = 0; i < need2load.size(); ++i) {
            MITK_INFO << "Loading " + need2load.at(i) + " files";

            if (need2load.at(i).compare("Pre", Qt::CaseSensitive)==0)
                prodPath = ScarCalculationsView::predir + mitk::IOUtil::GetDirectorySeparator();
            else
                prodPath = ScarCalculationsView::postdir + mitk::IOUtil::GetDirectorySeparator();

            ifstream prodFileRead;
            ofstream prodFileWrite;
            prodFileRead.open((prodPath + "prodThresholds.txt").toStdString());
            prodFileWrite.open((prodPathOut + "prodThresholds"+need2load.at(i)+".txt").toStdString());

            MITK_INFO << "READ FILE: " + prodPath + "prodThresholds.txt";
            MITK_INFO << "WRITE FILE: " + prodPathOut + "prodThresholds"+need2load.at(i)+".txt";

            for(int i = 0; i < 5; i++) {
                prodFileRead >> data1[i];
                prodFileWrite << data1[i] << "\n";
            }
            prodFileRead.close();
            prodFileWrite.close();
        }
    }
}

void ScarCalculationsView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.fandi_t1, SIGNAL(clicked()), this, SLOT(DoImageProcessing()));
    connect(m_Controls.fandi_t2, SIGNAL(clicked()), this, SLOT(GapMeasurement()));
    connect(m_Controls.fandi_t2_visualise, SIGNAL(currentIndexChanged(const QString&)), this, SLOT(GapMeasurementVisualisation(const QString&)));
    connect(m_Controls.fandi_t3, SIGNAL(clicked()), this, SLOT(BeforeAndAfterComp()));
    connect(m_Controls.fandi_t3_visualise, SIGNAL(clicked()), this, SLOT(BeforeAndAfterCompVisualisation()));
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(Sphericity()));
    connect(m_Controls.button_thres, SIGNAL(clicked()), this, SLOT(EditThreshold()));
    connect(m_Controls.button_saveth, SIGNAL(clicked()), this, SLOT(SaveNewThreshold()));
    connect(m_Controls.button_cancel, SIGNAL(clicked()), this, SLOT(CancelThresholdEdit()));
    connect(m_Controls.combo_thres, SIGNAL(currentIndexChanged(const QString&)),
            this, SLOT(SetNewThreshold(const QString&)));
    connect(m_Controls.comboBox, SIGNAL(currentIndexChanged(const QString&)),
            this, SLOT(CtrlPrePostSelection(const QString&)));

    //Setup renderer
    surfActor = vtkSmartPointer<vtkActor>::New();
    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0,0,0);

    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow =
            vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    m_Controls.widget_1->SetRenderWindow(renderWindow);
    m_Controls.widget_1->GetRenderWindow()->AddRenderer(renderer);

    //Setup keyboard interactor
    callBack = vtkSmartPointer<vtkCallbackCommand>::New();
    callBack->SetCallback(KeyCallBackFunc);
    callBack->SetClientData(this);
    interactor = m_Controls.widget_1->GetRenderWindow()->GetInteractor();
    interactor->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());
    interactor->GetInteractorStyle()->KeyPressActivationOff();
    interactor->GetInteractorStyle()->AddObserver(vtkCommand::KeyPressEvent, callBack);

    // Initialise
    iniPreSurf();

    if (surface.IsNotNull()) {
        QString prodPathOut = ScarCalculationsView::advdir + mitk::IOUtil::GetDirectorySeparator();
        pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
        pickedSeedIds->Initialize();
        pickedLineSeeds = vtkSmartPointer<vtkPolyData>::New();
        pickedLineSeeds->Initialize();
        pickedLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
        pickedCutterSeeds = vtkSmartPointer<vtkPolyData>::New();
        pickedCutterSeeds->Initialize();
        pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());

        csadv = std::unique_ptr<CemrgScarAdvanced>(new CemrgScarAdvanced());
        outprefix = "pre";

        csadv->SetOutputFileName((prodPathOut+outprefix+"encirclement.csv").toStdString());
        csadv->SetOutputPath(prodPathOut.toStdString());

        csadv->SetSurfaceAreaFilename("SurfaceAreaResults.txt");
        csadv->SetGapsFilename("GapsMeasurementsResults.txt");
        csadv->SetComparisonFilename("ComparisonResults.txt");

        csadv->SetOutputPrefix(outprefix.toStdString());
        csadv->SetFillThreshold(thres);
        csadv->SetInputData(surface->GetVtkPolyData());

        Visualiser();
        csadv->SetMaxScalar(this->maxScalar);
        SetShortcutLegend();
    }

    m_Controls.fandi_t2_visualise->addItem("");
    m_Controls.fandi_t2_visualise->setEnabled(false);
    m_Controls.fandi_t2_visualise->setVisible(false);
    m_Controls.fandi_t3_visualise->setEnabled(false);
    m_Controls.fandi_t3_visualise->setVisible(false);
    m_Controls.comboBox->setEnabled(true);
    m_Controls.button_1->setVisible(false);
    m_Controls.combo_thres->setEnabled(false);
    m_Controls.button_saveth->setEnabled(false);
    m_Controls.button_cancel->setEnabled(false);
}

void ScarCalculationsView::SetFocus() {

    m_Controls.fandi_t1->setFocus();
}

void ScarCalculationsView::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void ScarCalculationsView::iniPreSurf() {

    // Check for folders to exist!
    QDir pred(ScarCalculationsView::predir);
    QDir postd(ScarCalculationsView::postdir);
    QDir advd(ScarCalculationsView::advdir);

    if (!pred.exists() || !postd.exists()) {
        QMessageBox::warning(NULL, "Attention - Check folders' names!",
                             "The patient's folder does not seem to have the correct structure!");
        MITK_WARN << "The patient's folder does not seem to have the correct structure!";
        this->GetSite()->GetPage()->ResetPerspective();
        return;
    } else if (advd.mkdir(ScarCalculationsView::advdir)) {
        QMessageBox::warning(NULL, "Advanced analysis folder",
                             ("The advanced analysis folder was created in:\n\n"+
                              ScarCalculationsView::advdir.toStdString()).c_str());
        MITK_INFO << "The advanced analysis folder was created.";
    }
    // Load file information
    MITK_INFO << "Loading threshold information from file";
    double datainfo[5];
    ifstream prodFileRead;
    QString prodPathAdv = ScarCalculationsView::advdir + mitk::IOUtil::GetDirectorySeparator();
    prodFileRead.open((prodPathAdv + "prodThresholdsPre.txt").toStdString());

    MITK_INFO << "Read file: " + (prodPathAdv + "prodThresholdsPre.txt").toStdString();

    for(int i = 0; i < 5; i++) {
        prodFileRead >> datainfo[i];
    }
    prodFileRead.close();
    value = datainfo[0];
    method = datainfo[1];
    mean = datainfo[2];
    stdv = datainfo[3];
    thres = datainfo[4];

    // Convert to point data
    QString prename = ScarCalculationsView::preScarFile.isEmpty() ? "MaxScar.vtk" : ScarCalculationsView::preScarFile;
    QString shellPathPre = ScarCalculationsView::predir +  mitk::IOUtil::GetDirectorySeparator() + prename;
    MITK_INFO << "Shell PRE: " + shellPathPre.toStdString();
    mitk::Surface::Pointer shellpre = mitk::IOUtil::Load<mitk::Surface>(shellPathPre.toStdString());
    vtkSmartPointer<vtkCellDataToPointData> cell_to_point = vtkSmartPointer<vtkCellDataToPointData>::New();
    cell_to_point->SetInputData(shellpre->GetVtkPolyData());
    cell_to_point->PassCellDataOn();
    cell_to_point->Update();
    shellpre->SetVtkPolyData(cell_to_point->GetPolyDataOutput());
    mitk::IOUtil::Save(shellpre, (prodPathAdv+"MaxScarPre.vtk").toStdString());

    QString postname = ScarCalculationsView::postScarFile.isEmpty() ? "MaxScar.vtk" : ScarCalculationsView::postScarFile;
    QString shellPathPost = ScarCalculationsView::postdir +  mitk::IOUtil::GetDirectorySeparator() + postname;
    MITK_INFO << "Shell POST: " + shellPathPost.toStdString();
    mitk::Surface::Pointer shellpost = mitk::IOUtil::Load<mitk::Surface>(shellPathPost.toStdString());
    vtkSmartPointer<vtkCellDataToPointData> cell_to_point2 = vtkSmartPointer<vtkCellDataToPointData>::New();
    cell_to_point2->SetInputData(shellpost->GetVtkPolyData());
    cell_to_point2->PassCellDataOn();
    cell_to_point2->Update();
    shellpost->SetVtkPolyData(cell_to_point2->GetPolyDataOutput());
    mitk::IOUtil::Save(shellpost, (prodPathAdv+"MaxScarPost.vtk").toStdString());

    // Load preablation
    surface = shellpre;
    outprefix = "pre";

    QFileInfo txMaxScarPost(prodPathAdv+"MaxScarPost_Aligned.vtk");
    if (!txMaxScarPost.exists()) {
        MITK_INFO << "Transformed POST-ablation shell (MaxScarPost_Aligned.vtk) not found. Creating...";
        this->TransformMeshesForComparison();
    }
    m_Controls.comboBox->addItem("POST (ALIGNED)");
}

void ScarCalculationsView::KeyCallBackFunc(
        vtkObject*, long unsigned int, void* ClientData, void*) {

    ScarCalculationsView* self;
    self = reinterpret_cast<ScarCalculationsView*>(ClientData);
    std::string key = self->interactor->GetKeySym();
    vtkSmartPointer<vtkPolyData> poly_data = vtkSmartPointer<vtkPolyData>::New();
    poly_data = self->surface->GetVtkPolyData();

    if (key == "space") {
        MITK_INFO << "[INFO][KeyCallBackFunc] Pressed SPACE key.\n";
        self->PickCallBack();

    } else if (key == "Delete") {
        MITK_INFO << "[INFO][KeyCallBackFunc] Pressed DELETE key.\n";
        vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkPoints> points = self->pickedLineSeeds->GetPoints();
        for (int i=0; i<points->GetNumberOfPoints()-1; i++)
            newPoints->InsertNextPoint(points->GetPoint(i));
        self->pickedLineSeeds->SetPoints(newPoints);
        if (self->pickedSeedLabels.empty() == false)
            self->pickedSeedLabels.pop_back();
        self->m_Controls.widget_1->GetRenderWindow()->Render();

    } else if (key == "r" || key == "R") {
        MITK_INFO << "[INFO][KeyCallBackFunc] Pressed R key.\n";
        //Clear renderer
        self->renderer->RemoveAllViewProps();
        self->dijkstraActors.clear();
        self->pickedSeedLabels.clear();
        self->pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
        self->pickedSeedIds->Initialize();
        self->pickedLineSeeds = vtkSmartPointer<vtkPolyData>::New();
        self->pickedLineSeeds->Initialize();
        self->pickedLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
        self->pickedCutterSeeds = vtkSmartPointer<vtkPolyData>::New();
        self->pickedCutterSeeds->Initialize();
        self->pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());

        self->csadv->ResetValues();
        self->Visualiser();
        self->csadv->SetMaxScalar(self->maxScalar);
        self->SetShortcutLegend();
    }
}

void ScarCalculationsView::BinVisualiser() {

    MITK_INFO << "Binary Visualiser";
    double max_scalar=-1, min_scalar=1e9, s;
    int numlabels=2;
    vtkIntArray *scalars = vtkIntArray::New();
    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
    scalars = vtkIntArray::SafeDownCast(surface->GetVtkPolyData()->GetPointData()->GetScalars());

    for (vtkIdType i=0;i<surface->GetVtkPolyData()->GetNumberOfPoints();i++) {
        s = scalars->GetTuple1(i);
        if (s > max_scalar)
            max_scalar = s;
        if (s < min_scalar)
            min_scalar = s;
    }
    if (max_scalar==3){
        numlabels = 4;
    }

    vtkSmartPointer<vtkPolyDataMapper> surfMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    surfMapper->SetInputData(surface->GetVtkPolyData());
    surfMapper->SetScalarRange(min_scalar, max_scalar);
    surfMapper->SetScalarModeToUsePointData();
    surfMapper->ScalarVisibilityOn();
    lut->SetNumberOfTableValues(numlabels);
    lut->SetTableRange(min_scalar, max_scalar);
    lut->SetHueRange(0.6, 0.0);  // this is the way_neighbourhood_size you tell which colors you want to be displayed.
    lut->Build();     // this is important

    vtkSmartPointer<vtkScalarBarActor> scalarBar =
            vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable(surfMapper->GetLookupTable());
    scalarBar->SetWidth(0.75);
    scalarBar->SetHeight(0.3);
    scalarBar->SetTextPositionToPrecedeScalarBar();
    scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    scalarBar->GetPositionCoordinate()->SetValue( 0.9, 0.01 );
    scalarBar->SetNumberOfLabels(numlabels);

    surfMapper->SetLookupTable(lut);
    scalarBar->SetLookupTable(lut);

    vtkSmartPointer<vtkActor> surfActor = vtkSmartPointer<vtkActor>::New();
    surfActor->SetMapper(surfMapper);
    surfActor->GetProperty()->SetOpacity(1);

    renderer->AddActor(surfActor);
    renderer->AddActor2D(scalarBar);
}

void ScarCalculationsView::Visualiser() {

    MITK_INFO << "Visualiser";
    double max_scalar=-1, min_scalar=1e9,s;
    vtkFloatArray *scalars = vtkFloatArray::New();
    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
    scalars = vtkFloatArray::SafeDownCast(surface->GetVtkPolyData()->GetPointData()->GetScalars());
    for (vtkIdType i=0;i<surface->GetVtkPolyData()->GetNumberOfPoints();i++) {
        s = scalars->GetTuple1(i);
        if (s > max_scalar)
            max_scalar = s;
        if (s < min_scalar)
            min_scalar = s;
    }
    this->maxScalar = max_scalar;
    this->minScalar = min_scalar;

    vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    vtkSmartPointer<vtkSphereSource> glyphSource = vtkSmartPointer<vtkSphereSource>::New();
    glyph3D->SetInputData(pickedLineSeeds);
    glyph3D->SetSourceConnection(glyphSource->GetOutputPort());
    glyph3D->SetScaleModeToDataScalingOff();
    glyph3D->SetScaleFactor(surface->GetVtkPolyData()->GetLength()*0.01);
    glyph3D->Update();

    //Create a mapper and actor for glyph
    vtkSmartPointer<vtkPolyDataMapper> glyphMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    glyphMapper->SetInputConnection(glyph3D->GetOutputPort());
    vtkSmartPointer<vtkActor> glyphActor = vtkSmartPointer<vtkActor>::New();
    glyphActor->SetMapper(glyphMapper);
    glyphActor->GetProperty()->SetColor(1.0,0.0,0.0);
    glyphActor->PickableOff();
    renderer->AddActor(glyphActor);

    //Create a mapper and actor for surface
    vtkSmartPointer<vtkPolyDataMapper> surfMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    surfMapper->SetInputData(surface->GetVtkPolyData());
    surfMapper->SetScalarRange(min_scalar, max_scalar);
    surfMapper->SetScalarModeToUsePointData();
    surfMapper->ScalarVisibilityOn();
    lut->SetTableRange(min_scalar, max_scalar);
    lut->SetHueRange(0.3, 0.0);  // this is the way_neighbourhood_size you tell which colors you want to be displayed.
    lut->Build();     // this is important

    vtkSmartPointer<vtkScalarBarActor> scalarBar =
            vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable(surfMapper->GetLookupTable());
    scalarBar->SetWidth(0.075);
    scalarBar->SetHeight(0.3);
    scalarBar->SetTextPositionToPrecedeScalarBar();
    scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    scalarBar->GetPositionCoordinate()->SetValue( 0.9, 0.01 );
    scalarBar->SetNumberOfLabels(5);

    std::string scalarBarTitle = "Raw Signal \n Intensity";
    scalarBar->SetTitle(scalarBarTitle.c_str());
    scalarBar->SetVerticalTitleSeparation(15);

    surfMapper->SetLookupTable(lut);
    scalarBar->SetLookupTable(lut);

    vtkSmartPointer<vtkActor> surfActor = vtkSmartPointer<vtkActor>::New();
    surfActor->SetMapper(surfMapper);
    surfActor->GetProperty()->SetOpacity(1);
    renderer->AddActor(surfActor);
    renderer->AddActor2D(scalarBar);
}

void ScarCalculationsView::CtrlPrePostSelection(const QString& text) {

    QString cb = text;//m_Controls.comboBox->currentText();
    QMessageBox::warning(NULL, "Attention",
                         "Changing to: " + cb + "-ablation data");
    // Load file information
    MITK_INFO << "Loading threshold information from file";
    double datainfo[5];
    ifstream prodFileRead;
    QString prodPathAdv = ScarCalculationsView::advdir + mitk::IOUtil::GetDirectorySeparator();
    QString shellpath = prodPathAdv;

    if (cb.contains("PRE", Qt::CaseInsensitive)) {
        outprefix = "pre";
        if (cb.contains("TRANSFORMED", Qt::CaseInsensitive)) {
            shellpath = shellpath + "MaxScarPre_OnPost.vtk";
            outprefix = "Tx_" + outprefix;
        } else
            shellpath = shellpath + "MaxScarPre.vtk";
        MITK_INFO << "Changing to PRE-ablation scar map." + shellpath.toStdString();
        prodFileRead.open((prodPathAdv + "prodThresholdsPre.txt").toStdString());
    } else { //POST
        outprefix = "post";
        if (cb.contains("ALIGN", Qt::CaseInsensitive)) {
            shellpath = shellpath + "MaxScarPost_Aligned.vtk";
            outprefix = "Tx_" + outprefix;
        } else
            shellpath = shellpath + "MaxScarPost.vtk";
        MITK_INFO << "Changing to POST-ablation scar map." + shellpath.toStdString();
        prodFileRead.open((prodPathAdv + "prodThresholdsPost.txt").toStdString());
    }

    csadv->SetOutputFileName((prodPathAdv+outprefix+"encirclement.csv").toStdString());
    csadv->SetOutputPrefix(outprefix.toStdString());

    for(int i = 0; i < 5; i++) {
        prodFileRead >> datainfo[i];
    }
    prodFileRead.close();
    value = datainfo[0];
    method = datainfo[1];
    mean = datainfo[2];
    stdv = datainfo[3];
    thres = datainfo[4];

    mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(shellpath.toStdString());
    surface = shell;

    //Clear renderer
    renderer->RemoveAllViewProps();
    dijkstraActors.clear();
    pickedSeedLabels.clear();
    pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
    pickedSeedIds->Initialize();
    pickedLineSeeds = vtkSmartPointer<vtkPolyData>::New();
    pickedLineSeeds->Initialize();
    pickedLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
    pickedCutterSeeds = vtkSmartPointer<vtkPolyData>::New();
    pickedCutterSeeds->Initialize();
    pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());

    csadv->ResetValues();
    csadv->SetFillThreshold(thres);
    csadv->SetInputData(surface->GetVtkPolyData());
    MITK_INFO << "Number of points in shell: ";
    MITK_INFO << surface->GetVtkPolyData()->GetNumberOfPoints();

    Visualiser();
    csadv->SetMaxScalar(this->maxScalar);
    m_Controls.widget_1->GetRenderWindow()->Render();
    SetShortcutLegend();
}

void ScarCalculationsView::EditThreshold() {

    if (m_Controls.fandi_t2_visualise->isEnabled() || m_Controls.fandi_t3_visualise->isEnabled()) {
        QMessageBox::warning(NULL, "ATTENTION", "Press the Cancel button first.");
        return;
    }

    MITK_INFO << "Edit threshold button...";

    if (m_Controls.combo_thres->count()==0) { // sanity checks
        if (method==2)
            m_Controls.combo_thres->addItems({"1", "2", "2.3", "3.3", "4", "5"});
        else
            m_Controls.combo_thres->addItems({"0.86","0.97", "1.16", "1.2", "1.32"});
    }

    m_Controls.combo_thres->setEnabled(true);
    m_Controls.button_saveth->setEnabled(true);
    m_Controls.button_cancel->setEnabled(true);
    m_Controls.comboBox->setEnabled(false);

    std::string threShellPath = csadv->ThresholdedShell(thres);
    MITK_INFO << threShellPath;

    mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(threShellPath);
    surface = shell;

    //Clear renderer
    renderer->RemoveAllViewProps();
    dijkstraActors.clear();
    pickedSeedLabels.clear();
    pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
    pickedSeedIds->Initialize();
    pickedLineSeeds = vtkSmartPointer<vtkPolyData>::New();
    pickedLineSeeds->Initialize();
    pickedLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
    pickedCutterSeeds = vtkSmartPointer<vtkPolyData>::New();
    pickedCutterSeeds->Initialize();
    pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());

    csadv->ResetValues();
    BinVisualiser();
    m_Controls.widget_1->GetRenderWindow()->Render();
}

void ScarCalculationsView::SetNewThreshold(const QString& text) {

    if (m_Controls.combo_thres->count()==0)
        return;
    MITK_INFO << "New threshold selected...";

    QString prodPath = ScarCalculationsView::advdir + mitk::IOUtil::GetDirectorySeparator();
    QString cb = m_Controls.comboBox->currentText();
    QString outname = "prodThresholds";

    if (cb.contains("PRE", Qt::CaseInsensitive))
        outname = outname + "Pre.txt";
    else//POST
        outname = outname + "Post.txt";

    MITK_INFO << "Reading thresholds from: " + prodPath + outname;

    value = text.toDouble();
    thres = (method == 1) ? mean*value : mean+value*stdv;

    std::string threShellPath = csadv->ThresholdedShell(thres);

    mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(threShellPath);
    surface = shell;

    //Clear renderer
    renderer->RemoveAllViewProps();
    dijkstraActors.clear();
    pickedSeedLabels.clear();
    pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
    pickedSeedIds->Initialize();
    pickedLineSeeds = vtkSmartPointer<vtkPolyData>::New();
    pickedLineSeeds->Initialize();
    pickedLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
    pickedCutterSeeds = vtkSmartPointer<vtkPolyData>::New();
    pickedCutterSeeds->Initialize();
    pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());

    csadv->ResetValues();
    BinVisualiser();
    m_Controls.widget_1->GetRenderWindow()->Render();
}

void ScarCalculationsView::SaveNewThreshold() {

    m_Controls.combo_thres->setEnabled(false);
    m_Controls.button_saveth->setEnabled(false);
    m_Controls.button_cancel->setEnabled(false);
    m_Controls.comboBox->setEnabled(true);

    QString prodPath = ScarCalculationsView::advdir + mitk::IOUtil::GetDirectorySeparator();
    QString shellpath = prodPath;
    QString cb = m_Controls.comboBox->currentText();
    QString outname = "prodThresholds";

    if (cb.contains("PRE", Qt::CaseInsensitive)) {
        outname = outname + "Pre.txt";
        outprefix = "pre";
        if (cb.contains("TRANSFORMED", Qt::CaseInsensitive)) {
            shellpath = prodPath + "MaxScarPre_OnPost.vtk";
            outprefix = "Tx_" + outprefix;
        } else
            shellpath = prodPath + "MaxScarPre.vtk";
    }
    else {//POST
        outname = outname + "Post.txt";
        outprefix = "post";
        if (cb.contains("ALIGN", Qt::CaseInsensitive)) {
            shellpath = prodPath + "MaxScarPost_Aligned.vtk";
            outprefix = "Tx_" + outprefix;
        } else
            shellpath = prodPath + "MaxScarPost.vtk";
    }
    csadv->SetOutputPrefix(outprefix.toStdString());

    MITK_INFO << "Writing threshold information to: " + prodPath + outname;
    SetThresholdValuesToFile(prodPath + outname);
    GetThresholdValuesFromFile(prodPath + outname);

    mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(shellpath.toStdString());
    surface = shell;

    //Clear renderer
    renderer->RemoveAllViewProps();
    dijkstraActors.clear();
    pickedSeedLabels.clear();
    pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
    pickedSeedIds->Initialize();
    pickedLineSeeds = vtkSmartPointer<vtkPolyData>::New();
    pickedLineSeeds->Initialize();
    pickedLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
    pickedCutterSeeds = vtkSmartPointer<vtkPolyData>::New();
    pickedCutterSeeds->Initialize();
    pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());

    csadv->ResetValues();
    csadv->SetFillThreshold(thres);
    csadv->SetInputData(surface->GetVtkPolyData());

    Visualiser();
    csadv->SetMaxScalar(this->maxScalar);
    m_Controls.widget_1->GetRenderWindow()->Render();
    SetShortcutLegend();
}

void ScarCalculationsView::CancelThresholdEdit() {
    m_Controls.combo_thres->setEnabled(false);
    m_Controls.button_saveth->setEnabled(false);
    m_Controls.button_cancel->setEnabled(false);
    m_Controls.comboBox->setEnabled(true);

    if (m_Controls.fandi_t2_visualise->isEnabled()) {
        m_Controls.fandi_t2_visualise->setEnabled(false);
        m_Controls.fandi_t2_visualise->setVisible(false);
    }
    if (m_Controls.fandi_t3_visualise->isEnabled()) {
        m_Controls.fandi_t3_visualise->setEnabled(false);
        m_Controls.fandi_t3_visualise->setVisible(false);
    }

    QString prodPath = ScarCalculationsView::advdir + mitk::IOUtil::GetDirectorySeparator();
    QString shellpath = prodPath;
    QString cb = m_Controls.comboBox->currentText();
    QString outname = "prodThresholds";

    if (cb.contains("PRE", Qt::CaseInsensitive)) {
        outname = outname + "Pre.txt";
        outprefix = "pre";
        if (cb.contains("TRANSFORMED", Qt::CaseInsensitive)) {
            shellpath = prodPath + "MaxScarPre_OnPost.vtk";
            outprefix = "Tx_" + outprefix;
        } else
            shellpath = prodPath + "MaxScarPre.vtk";
    }
    else {//POST
        outname = outname + "Post.txt";
        outprefix = "post";
        if (cb.contains("ALIGN", Qt::CaseInsensitive)) {
            shellpath = prodPath + "MaxScarPost_Aligned.vtk";
            outprefix = "Tx_" + outprefix;
        } else
            shellpath = prodPath + "MaxScarPost.vtk";
    }
    csadv->SetOutputPrefix(outprefix.toStdString());

    GetThresholdValuesFromFile(prodPath + outname);

    mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(shellpath.toStdString());
    surface = shell;

    //Clear renderer
    renderer->RemoveAllViewProps();
    dijkstraActors.clear();
    pickedSeedLabels.clear();
    pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
    pickedSeedIds->Initialize();
    pickedLineSeeds = vtkSmartPointer<vtkPolyData>::New();
    pickedLineSeeds->Initialize();
    pickedLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
    pickedCutterSeeds = vtkSmartPointer<vtkPolyData>::New();
    pickedCutterSeeds->Initialize();
    pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());

    csadv->ResetValues();
    csadv->SetFillThreshold(thres);
    csadv->SetInputData(surface->GetVtkPolyData());

    Visualiser();
    csadv->SetMaxScalar(this->maxScalar);
    m_Controls.widget_1->GetRenderWindow()->Render();
    SetShortcutLegend();
}

void ScarCalculationsView::PickCallBack() {

    vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(1E-4 * surface->GetVtkPolyData()->GetLength());
    int* eventPosition = interactor->GetEventPosition();
    int result = picker->Pick(float(eventPosition[0]), float(eventPosition[1]), 0.0, renderer);
    if (result == 0) return;
    double* pickPosition = picker->GetPickPosition();
    vtkIdList* pickedCellPointIds = surface->GetVtkPolyData()->GetCell(picker->GetCellId())->GetPointIds();

    double distance;
    int pickedSeedId = -1;
    double minDistance = 1E10;
    for (int i=0; i<pickedCellPointIds->GetNumberOfIds(); i++) {
        distance = vtkMath::Distance2BetweenPoints(
                    pickPosition, surface->GetVtkPolyData()->GetPoint(pickedCellPointIds->GetId(i)));
        if (distance < minDistance) {
            minDistance = distance;
            pickedSeedId = pickedCellPointIds->GetId(i);
        }//_if
    }//_for
    if (pickedSeedId == -1)
        pickedSeedId = pickedCellPointIds->GetId(0);

    pickedSeedIds->InsertNextId(pickedSeedId);
    double* point = surface->GetVtkPolyData()->GetPoint(pickedSeedId);
    pickedLineSeeds->GetPoints()->InsertNextPoint(point);
    pickedLineSeeds->Modified();
    m_Controls.widget_1->GetRenderWindow()->Render();
}

// Button functionalities
void ScarCalculationsView::DoImageProcessing() {

    if (m_Controls.fandi_t2_visualise->isEnabled() || m_Controls.fandi_t3_visualise->isEnabled()) {
        QMessageBox::warning(NULL, "ATTENTION", "Press the Cancel button first.");
        return;
    }
    if (m_Controls.button_saveth->isEnabled()) {
        QMessageBox::warning(NULL, "ATTENTION", "Press the Cancel button or save your threshold first.");
        return;
    }

    // F&I T1
    MITK_INFO << "[INFO] F&I Task 1. Comparison of Surface Area.";
    csadv->GetSurfaceAreaFromThreshold(thres, maxScalar);
    csadv->ScarScore(thres);
    QMessageBox::warning(NULL, "F&I T1 - FINISHED CALCULATION",
                         (csadv->PrintThresholdResults(mean, stdv, value)).c_str());
}

void ScarCalculationsView::GapMeasurement() {

    if (m_Controls.fandi_t2_visualise->isEnabled() || m_Controls.fandi_t3_visualise->isEnabled()) {
        QMessageBox::warning(NULL, "ATTENTION", "Press the Cancel button first.");
        return;
    }
    if (m_Controls.button_saveth->isEnabled()) {
        QMessageBox::warning(NULL, "ATTENTION", "Press the Cancel button or save your threshold first.");
        return;
    }

    // F&I T2
    if (pickedSeedIds->GetNumberOfIds()==0) {

        QMessageBox::warning(NULL, "Attention - No points selected.",
                             "Please select at least five points surrounding a vein.");
        MITK_WARN << "Please select at least five points surrounding a vein.";

    } else {

        MITK_INFO << "[F&I Task 2]. Measurement of ablation gaps.";
        MITK_INFO << "Passing selected IDs to underlying functionalities.";
        int lim = this->pickedSeedIds->GetNumberOfIds();
        std::vector<int> v;
        for(int i=0; i<lim; i++){
            v.push_back(this->pickedSeedIds->GetId(i));
        }

        MITK_INFO << "Creating shortest path and corridor.";
        // Choose parameters: neighbourhood size, left/right prefix
        QDialog* inputs = new QDialog(0,0);
        m_UICorridor.setupUi(inputs);
        connect(m_UICorridor.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
        connect(m_UICorridor.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
        int dialogCode = inputs->exec();

        //Act on dialog return code
        if (dialogCode == QDialog::Accepted) {

            bool ok1;
            int thickness = m_UICorridor.txtbox_thick->text().toInt(&ok1);
            std::string lrpre = m_UICorridor.rb_left->isChecked() ? "left_" : "right_";

            csadv->SetWeightedCorridorBool(m_UICorridor.checkBox_weighted->isChecked());

            if (!ok1) { //Set default values
                QMessageBox::warning(NULL, "Attention", "Using default thickness!");
                thickness = 3;
            }
            // Calculate neighbourhood size from thickness
            MITK_INFO << "Neighbourhood size:" + csadv->num2str(thickness,0);
            MITK_INFO << "Side: " + lrpre;
            csadv->SetNeighbourhoodSize(thickness);
            csadv->SetLeftRightPrefix(lrpre);
            csadv->CorridorFromPointList(v);
            if (m_Controls.fandi_t2_visualise->findText(QString::fromStdString(csadv->GetPrefix())+"exploration_scalars")==-1)
                m_Controls.fandi_t2_visualise->addItem(QString::fromStdString(csadv->GetPrefix())+"exploration_scalars");
            csadv->ClearLeftRightPrefix();
            this->dijkstraActors = csadv->GetPathsMappersAndActors();

            for (int i=0;(unsigned)i<this->dijkstraActors.size();i++){
                renderer->AddActor(this->dijkstraActors[i]);
            }

            m_Controls.widget_1->GetRenderWindow()->Render();

            QMessageBox::warning(NULL, "F&I T2 - FINISHED CALCULATION",
                                 (csadv->PrintAblationGapsResults(mean, stdv, value)).c_str());
            csadv->ResetValues();
            inputs->deleteLater();

        } else if (dialogCode == QDialog::Rejected) {
            QMessageBox::warning(NULL, "F&I T2 - CALCULATION CANCELLED",
                                 "'Cancel' button pressed, no calculations were made.");
            inputs->close();
            inputs->deleteLater();
        }//_if
    }

    m_Controls.fandi_t2_visualise->setEnabled(true);
    m_Controls.fandi_t2_visualise->setVisible(true);
    m_Controls.button_cancel->setEnabled(true);
}

void ScarCalculationsView::GapMeasurementVisualisation(const QString& text) {

    MITK_INFO << "Visualisation of exploiration corridor";
    QString cb = text;

    if (!cb.isEmpty()) {

        MITK_INFO << ("Current text: " + cb).toStdString();
        QString prodPath = ScarCalculationsView::advdir + mitk::IOUtil::GetDirectorySeparator();
        QFileInfo fi(prodPath + cb + ".vtk");
        MITK_INFO << ("Changing to file" + fi.absoluteFilePath()).toStdString();

        if (fi.exists()) {

            mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(fi.absoluteFilePath().toStdString());
            surface = shell;

            //Clear renderer
            renderer->RemoveAllViewProps();
            dijkstraActors.clear();
            pickedSeedLabels.clear();
            pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
            pickedSeedIds->Initialize();
            pickedLineSeeds = vtkSmartPointer<vtkPolyData>::New();
            pickedLineSeeds->Initialize();
            pickedLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
            pickedCutterSeeds = vtkSmartPointer<vtkPolyData>::New();
            pickedCutterSeeds->Initialize();
            pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());

            csadv->ResetValues();
            BinVisualiser();
            m_Controls.widget_1->GetRenderWindow()->Render();
        }
    } else
        MITK_INFO << "Empty name of file. ";
}

void ScarCalculationsView::BeforeAndAfterComp() {

    if (m_Controls.fandi_t2_visualise->isEnabled() || m_Controls.fandi_t3_visualise->isEnabled()) {
        QMessageBox::warning(NULL, "ATTENTION", "Press the Cancel button first.");
        return;
    }
    if (m_Controls.button_saveth->isEnabled()) {
        QMessageBox::warning(NULL, "ATTENTION", "Press the Cancel button or save your threshold first.");
        return;
    }

    // F&I T3
    MITK_INFO << "[INFO] F&I Task 3. Measurement of scar overlap.\n";
    QString current = m_Controls.comboBox->currentText();
    double valpre, valpost;

    QString outpath = ScarCalculationsView::advdir + mitk::IOUtil::GetDirectorySeparator();
    QString outScarMap = outpath + "MaxScarPost_Aligned.vtk";

    QFileInfo txMaxScarPost(outScarMap);
    if (!txMaxScarPost.exists()) {
        MITK_INFO << "Transformed POST-ablation shell (MaxScarPost_Aligned.vtk) not found. Creating...";
        this->TransformMeshesForComparison();
    }

    QString preShellPath = outpath + "MaxScarPre.vtk";
    QString preThresPath = outpath + "prodThresholdsPre.txt";
    mitk::Surface::Pointer shellpre = mitk::IOUtil::Load<mitk::Surface>(preShellPath.toStdString());

    GetThresholdValuesFromFile(preThresPath);
    valpre = value;
    csadv->SetInputData(shellpre->GetVtkPolyData());
    csadv->SetOutputPrefix("pre");
    csadv->GetSurfaceAreaFromThreshold(thres, maxScalar);
    csadv->ScarScore(thres);

    QString postThresPath = outpath + "prodThresholdsPost.txt";
    mitk::Surface::Pointer shellpost = mitk::IOUtil::Load<mitk::Surface>(outScarMap.toStdString());

    GetThresholdValuesFromFile(postThresPath);
    valpost = value;
    csadv->SetInputData(shellpost->GetVtkPolyData());
    csadv->SetOutputPrefix("post");
    csadv->GetSurfaceAreaFromThreshold(thres, maxScalar);
    csadv->ScarScore(thres);

    this->CopyScalarValues();
    if (m_Controls.comboBox->findText("PRE (TRANSFORMED)", Qt::MatchExactly)==-1)
        m_Controls.comboBox->addItem("PRE (TRANSFORMED)");

    QMessageBox::warning(NULL, "F&I T3 - FINISHED CALCULATION",
                         (csadv->PrintScarOverlapResults(valpre, valpost)).c_str());
    // back to normal
    this->CtrlPrePostSelection(current);
    m_Controls.fandi_t3_visualise->setEnabled(true);
    m_Controls.fandi_t3_visualise->setVisible(true);
    m_Controls.button_cancel->setEnabled(true);
}

void ScarCalculationsView::BeforeAndAfterCompVisualisation() {

    MITK_INFO << "Visualisation of pre/post comparison.";
    QString outpath = ScarCalculationsView::advdir + mitk::IOUtil::GetDirectorySeparator();
    QString preMap = outpath + "MaxScarPre_OnPost.vtk";
    QString preThresPath = outpath + "prodThresholdsPre.txt";
    QString postMap = outpath + "MaxScarPost_Aligned.vtk";
    QString postThresPath = outpath + "prodThresholdsPost.txt";

    double prethresh, postthresh;

    mitk::Surface::Pointer presh = mitk::IOUtil::Load<mitk::Surface>(preMap.toStdString());
    mitk::Surface::Pointer postsh = mitk::IOUtil::Load<mitk::Surface>(postMap.toStdString());
    GetThresholdValuesFromFile(preThresPath);
    prethresh = thres;
    GetThresholdValuesFromFile(postThresPath);
    postthresh = thres;
    std::string overlapShellPath = csadv->ScarOverlap(presh->GetVtkPolyData(), prethresh, postsh->GetVtkPolyData(), postthresh);
    MITK_INFO << overlapShellPath;

    mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(overlapShellPath);
    surface = shell;

    //Clear renderer
    renderer->RemoveAllViewProps();
    dijkstraActors.clear();
    pickedSeedLabels.clear();
    pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
    pickedSeedIds->Initialize();
    pickedLineSeeds = vtkSmartPointer<vtkPolyData>::New();
    pickedLineSeeds->Initialize();
    pickedLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
    pickedCutterSeeds = vtkSmartPointer<vtkPolyData>::New();
    pickedCutterSeeds->Initialize();
    pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());

    csadv->ResetValues();
    BinVisualiser();
    m_Controls.widget_1->GetRenderWindow()->Render();
}

void ScarCalculationsView::Sphericity() {
}

// Helper functions
void ScarCalculationsView::TransformMeshesForComparison() {
    MITK_INFO << "[ATTENTION] Implementation of alignement routines.";

    bool successful;
    QString outpath = ScarCalculationsView::advdir + mitk::IOUtil::GetDirectorySeparator();
    QString sourcename = outpath + "MaxScarPost.vtk";
    QString targetname = outpath + "MaxScarPre.vtk";
    QString alignedname = outpath + "MaxScarPost_Aligned.vtk";
    QString testTX = outpath + "translation.dof";
    this->BusyCursorOn();
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->ExecuteSimpleTranslation(directory, sourcename, targetname, testTX);
    this->BusyCursorOff();
    if (cmd->IsOutputSuccessful(testTX)) {
        MITK_INFO << "[...] DOF file created successfully, attempting transformation.";
        this->BusyCursorOn();
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->ExecuteTransformationOnPoints(directory, sourcename, alignedname, testTX);
        this->BusyCursorOff();
        successful = cmd->IsOutputSuccessful(alignedname);
    } else {
        MITK_INFO << "[...] DOF file not created.";
        successful = false;
    }

    MITK_INFO << "Command Line Operations Finished!";
    MITK_WARN(!successful) << "Aligned file NOT created. Check the log.";
}

void ScarCalculationsView::GetThresholdValuesFromFile(QString filepath) {

    double data1[5];
    ifstream prodFileRead;
    prodFileRead.open(filepath.toStdString());

    MITK_INFO << "READ FILE: " + filepath;

    for(int i = 0; i < 5; i++)
        prodFileRead >> data1[i];

    value = data1[0];
    method = data1[1];
    mean = data1[2];
    stdv = data1[3];
    thres = data1[4];
    prodFileRead.close();
}

void ScarCalculationsView::SetThresholdValuesToFile(QString filepath) {

    double data1[5];
    ofstream prodFileWrite;
    prodFileWrite.open(filepath.toStdString());

    data1[0] = value;
    data1[1] = method;
    data1[2] = mean;
    data1[3] = stdv;
    data1[4] = thres;

    MITK_INFO << "WRITE FILE: " + filepath;

    for(int i = 0; i < 5; i++)
        prodFileWrite << data1[i] << "\n";
    prodFileWrite.close();
}

void ScarCalculationsView::SetShortcutLegend() {

    std::string mymethod = (method==1) ? "V*IIR" : "mean + V*stdv";
    mymethod = "METHOD: " + mymethod;

    std::string title = " "+ m_Controls.comboBox->currentText().toStdString() + "-Ablation ";
    title += "\n\n Value for threshold (V): " + csadv->num2str(value, 1);
    title += "\n "+ mymethod + " = " + csadv->num2str(thres, 1);
    title += "\n\n Shortcuts for gap measurement:\n ";

    std::string shortcuts = title + "Space: add seed point\n Delete: remove seed point\n R: Reset points.";

    vtkSmartPointer<vtkTextActor> txtActor = vtkSmartPointer<vtkTextActor>::New();
    txtActor->SetInput(shortcuts.c_str());
    txtActor->GetTextProperty()->SetFontSize(14);
    txtActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
    renderer->AddActor2D(txtActor);
}

void ScarCalculationsView::CopyScalarValues() {

    QString outpath = ScarCalculationsView::advdir + mitk::IOUtil::GetDirectorySeparator();
    QString sourcename = outpath + "MaxScarPost_Aligned.vtk";
    QString targetname = outpath + "MaxScarPre.vtk";

    mitk::Surface::Pointer _source = mitk::IOUtil::Load<mitk::Surface>(sourcename.toStdString());
    mitk::Surface::Pointer _target = mitk::IOUtil::Load<mitk::Surface>(targetname.toStdString());

    MITK_INFO << "[ATTENTION] Copying scalar values from MaxScarPre into MaxScarPost_Aligned";
    csadv->SetSourceAndTarget(_source->GetVtkPolyData(), _target->GetVtkPolyData());
    csadv->TransformSource2Target();
}
