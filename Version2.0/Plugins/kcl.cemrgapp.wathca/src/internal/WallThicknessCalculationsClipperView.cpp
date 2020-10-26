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
 * Morphological Quantification
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
#include <berryIWorkbenchWindow.h>
#include <berryIWorkbenchPage.h>

// Qmitk
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkNodePredicateProperty.h>
#include <mitkImage.h>
#include "WallThicknessCalculationsClipperView.h"
#include "WallThicknessCalculationsView.h"

// VTK
#include <vtkGlyph3D.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkProperty.h>
#include <vtkCellPicker.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCell.h>
#include <vtkMath.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkDecimatePro.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkPolyDataConnectivityFilter.h>

//ITK
#include <itkAddImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include "itkLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelSelectionLabelMapFilter.h"

// Qt
#include <QMessageBox>
#include <QDesktopWidget>

// CemrgAppModule
#include <CemrgAtriaClipper.h>
#include <CemrgCommandLine.h>

QString WallThicknessCalculationsClipperView::fileName;
QString WallThicknessCalculationsClipperView::directory;
const std::string WallThicknessCalculationsClipperView::VIEW_ID = "org.mitk.views.wathcaclipperview";

void WallThicknessCalculationsClipperView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(CtrLines()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(CtrPlanes()));
    connect(m_Controls.button_4, SIGNAL(clicked()), this, SLOT(ClipperImage()));
    connect(m_Controls.slider, SIGNAL(valueChanged(int)), this, SLOT(CtrPlanesPlacer()));
    connect(m_Controls.spinBox, SIGNAL(valueChanged(double)), this, SLOT(CtrPlanesPlacer()));
    connect(m_Controls.comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(CtrLinesSelector(int)));

    //Create GUI widgets
    inputs = new QDialog(0,0);
    m_Labels.setupUi(inputs);
    connect(m_Labels.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_Labels.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));

    //Setup renderer
    surfActor = vtkSmartPointer<vtkActor>::New();
    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0,0,0);
    vtkSmartPointer<vtkTextActor> txtActor = vtkSmartPointer<vtkTextActor>::New();
    std::string shortcuts = "R: reset centrelines\nSpace: add seed point\nDelete: remove seed point";
    txtActor->SetInput(shortcuts.c_str());
    txtActor->GetTextProperty()->SetFontSize(14);
    txtActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
    renderer->AddActor2D(txtActor);

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

    //Initialisation
    iniPreSurf();
    if (surface.IsNotNull()) {
        pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
        pickedSeedIds->Initialize();
        pickedLineSeeds = vtkSmartPointer<vtkPolyData>::New();
        pickedLineSeeds->Initialize();
        pickedLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
        pickedCutterSeeds = vtkSmartPointer<vtkPolyData>::New();
        pickedCutterSeeds->Initialize();
        pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
        clipper = std::unique_ptr<CemrgAtriaClipper>(new CemrgAtriaClipper(directory, surface));
        Visualiser();
    }
    m_Controls.slider->setEnabled(false);
    m_Controls.spinBox->setEnabled(false);
    m_Controls.comboBox->setEnabled(false);
}

void WallThicknessCalculationsClipperView::SetFocus() {

    m_Controls.button_1->setFocus();
}

void WallThicknessCalculationsClipperView::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*src*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

WallThicknessCalculationsClipperView::~WallThicknessCalculationsClipperView() {

    inputs->deleteLater();
}

void WallThicknessCalculationsClipperView::SetDirectoryFile(const QString directory, const QString fileName) {

    WallThicknessCalculationsClipperView::fileName = fileName;
    WallThicknessCalculationsClipperView::directory = directory;
}

void WallThicknessCalculationsClipperView::iniPreSurf() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation to clip!");
        this->GetSite()->GetPage()->ResetPerspective();
        return;
    }

    //Find the selected node
    QString path;
    mitk::DataNode::Pointer segNode = nodes.at(0);
    mitk::BaseData::Pointer data = segNode->GetData();
    if (data) {
        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {

            //Check seg node name
            if (segNode->GetName().compare(fileName.left(fileName.length()-4).toStdString()) != 0) {
                QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
                this->GetSite()->GetPage()->ResetPerspective();
                return;
            }//_if

            //Ask for user input to set the parameters
            QDialog* inputs = new QDialog(0,0);
            m_UIMeshing.setupUi(inputs);
            connect(m_UIMeshing.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
            connect(m_UIMeshing.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
            int dialogCode = inputs->exec();

            //Act on dialog return code
            if (dialogCode == QDialog::Accepted) {

                bool ok1, ok2, ok3, ok4, ok5;
                int iter = m_UIMeshing.lineEdit_1->text().toInt(&ok1);
                float th = m_UIMeshing.lineEdit_2->text().toFloat(&ok2);
                int blur = m_UIMeshing.lineEdit_3->text().toInt(&ok3);
                int smth = m_UIMeshing.lineEdit_4->text().toInt(&ok4);
                float ds = m_UIMeshing.lineEdit_5->text().toFloat(&ok5);

                //Set default values
                if (!ok1 || !ok2 || !ok3 || !ok4 || !ok5)
                    QMessageBox::warning(NULL, "Attention", "Reverting to default parameters!");
                if (!ok1) iter = 1;
                if (!ok2) th   = 0.5;
                if (!ok3) blur = 0;
                if (!ok4) smth = 0;
                if (!ok5) ds   = 0.99;
                //_if

                this->BusyCursorOn();
                path = directory + mitk::IOUtil::GetDirectorySeparator() + "temp.nii";
                mitk::IOUtil::Save(image, path.toStdString());
                mitk::ProgressBar::GetInstance()->AddStepsToDo(3);
                std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
                QString output = cmd->ExecuteSurf(directory, path, "close", iter, th, blur, smth);
                QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
                this->BusyCursorOff();

                //Decimate the mesh to visualise
                mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(output.toStdString());
                vtkSmartPointer<vtkDecimatePro> deci = vtkSmartPointer<vtkDecimatePro>::New();
                deci->SetInputData(shell->GetVtkPolyData());
                deci->SetTargetReduction(ds);
                deci->PreserveTopologyOn();
                deci->Update();
                vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
                connectivityFilter->SetInputConnection(deci->GetOutputPort());
                connectivityFilter->ColorRegionsOff();
                connectivityFilter->SetExtractionModeToLargestRegion();
                connectivityFilter->Update();
                vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
                smoother->SetInputConnection(connectivityFilter->GetOutputPort());
                smoother->SetNumberOfIterations(30);
                smoother->BoundarySmoothingOff();
                smoother->FeatureEdgeSmoothingOff();
                smoother->SetFeatureAngle(120.0);
                smoother->SetPassBand(.01);
                smoother->NonManifoldSmoothingOn();
                smoother->NormalizeCoordinatesOn();
                smoother->Update();
                vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
                cleaner->SetInputConnection(smoother->GetOutputPort());
                cleaner->PieceInvariantOn();
                cleaner->ConvertLinesToPointsOn();
                cleaner->ConvertStripsToPolysOn();
                cleaner->PointMergingOn();
                cleaner->Update();
                vtkSmartPointer<vtkPolyDataNormals> computeNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
                computeNormals->SetInputConnection(cleaner->GetOutputPort());
                computeNormals->SetFeatureAngle(360.0f);
                computeNormals->AutoOrientNormalsOn();
                computeNormals->FlipNormalsOff();
                computeNormals->Update();
                shell->SetVtkPolyData(computeNormals->GetOutput());
                surface = shell;

                //Tidy up data
                remove(path.toStdString().c_str());
                inputs->deleteLater();

            } else if (dialogCode == QDialog::Rejected) {
                inputs->close();
                inputs->deleteLater();
                this->GetSite()->GetPage()->ResetPerspective();
                return;
            }//_if

        } else {
            QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
            this->GetSite()->GetPage()->ResetPerspective();
            return;
        }//_image

    } else {
        this->GetSite()->GetPage()->ResetPerspective();
        return;
    }//_if
}

void WallThicknessCalculationsClipperView::CtrLines() {

    if (clipper->GetCentreLines().size() == 0 && pickedSeedLabels.size() == 0) {

        QMessageBox::warning(NULL, "Attention", "Please maske sure you have selected seeds!");
        return;

    } else {

        //Retrieve centrelines and labels
        if (clipper->GetCentreLines().size() != 0)
            QMessageBox::information(NULL, "Attention", "You are using precomputed centrelines!");

        this->BusyCursorOn();
        for (unsigned int i=0; i<clipper->GetCentreLines().size(); i++) {
            if (i == 0) pickedSeedLabels.clear();
            vtkSmartPointer<vtkPolyData> pd = clipper->GetCentreLines().at(i)->GetOutput();
            vtkSmartPointer<vtkIntArray> lb = vtkIntArray::SafeDownCast(pd->GetFieldData()->GetAbstractArray("PickedSeedLabels"));
            pickedSeedLabels.push_back(lb->GetValue(0));
        }//_for
        bool successful = clipper->ComputeCtrLines(pickedSeedLabels, pickedSeedIds, m_Controls.checkBox->isChecked());
        m_Controls.checkBox->setChecked(clipper->GetCentreLinesOrientation());
        this->BusyCursorOff();

        //Check for failure
        if (!successful) {
            QMessageBox::critical(NULL, "Attention", "Computation of Centrelines Failed!");
            return;
        }//_if
    }//_if

    //Set surface opacity
    renderer->RemoveAllViewProps();
    vtkSmartPointer<vtkPolyDataMapper> surfMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    surfMapper->SetInputData(surface->GetVtkPolyData());
    surfMapper->ScalarVisibilityOff();
    surfActor->SetMapper(surfMapper);
    surfActor->GetProperty()->SetOpacity(0.5);
    renderer->AddActor(surfActor);

    //Setup shortcut keys
    std::string sk;
    vtkSmartPointer<vtkTextActor> txtActor = vtkSmartPointer<vtkTextActor>::New();
    sk = "A: revert to automatic\nO: change opacity\nArrows: move planes\nSpace: add clip point\nDelete: remove clip point";
    txtActor->SetInput(sk.c_str());
    txtActor->GetTextProperty()->SetFontSize(14);
    txtActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
    renderer->AddActor2D(txtActor);

    //Create a mapper and actor for centre lines
    std::vector<vtkSmartPointer<vtkvmtkPolyDataCenterlines>> ctrLines = clipper->GetCentreLines();
    for (unsigned int i=0; i<ctrLines.size(); i++) {
        vtkSmartPointer<vtkPolyDataMapper> linesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        linesMapper->SetInputData(ctrLines.at(i)->GetOutput());
        linesMapper->ScalarVisibilityOff();
        vtkSmartPointer<vtkActor> linesActor = vtkSmartPointer<vtkActor>::New();
        linesActor->SetMapper(linesMapper);
        linesActor->GetProperty()->SetOpacity(1.0);
        linesActor->GetProperty()->SetColor(1,0,0);
        renderer->AddActor(linesActor);
    }//_for
    m_Controls.widget_1->GetRenderWindow()->Render();

    //Adjust controllers
    m_Controls.button_1->setEnabled(false);
    m_Controls.checkBox->setEnabled(false);
}

void WallThicknessCalculationsClipperView::CtrPlanes() {

    if (clipper->GetCentreLines().size() == 0 || pickedSeedLabels.size() == 0) {

        QMessageBox::warning(NULL, "Attention", "Please maske sure you have computed centre lines!");
        return;

    } else {

        this->BusyCursorOn();
        bool successful = clipper->ComputeCtrLinesClippers(pickedSeedLabels);
        this->BusyCursorOff();

        //Check for failure
        if (!successful) {
            QMessageBox::critical(NULL, "Attention", "Computation of Clipper Planes Failed!");
            return;
        }//_if
    }//_if

    std::vector<vtkSmartPointer<vtkRegularPolygonSource>> ctrPlanes = clipper->GetCentreLinePolyPlanes();
    std::vector<vtkSmartPointer<vtkvmtkPolyDataCenterlines>> ctrLines = clipper->GetCentreLines();
    vtkSmartPointer<vtkRegularPolygonSource> ctrPlane = ctrPlanes.at(0);
    vtkSmartPointer<vtkvmtkPolyDataCenterlines> ctrLine = ctrLines.at(0);
    int index = ctrLine->GetOutput()->FindPoint(ctrPlane->GetCenter());

    //Create a mapper and actor for centre lines clippers
    for (unsigned int i=0; i<ctrPlanes.size(); i++) {
        vtkSmartPointer<vtkPolyDataMapper> clipperMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        clipperMapper->SetInputConnection(ctrPlanes.at(i)->GetOutputPort());
        clipperMapper->ScalarVisibilityOff();
        vtkSmartPointer<vtkActor> clipperActor = vtkSmartPointer<vtkActor>::New();
        clipperActor->SetMapper(clipperMapper);
        clipperActor->GetProperty()->SetOpacity(1.0);
        clipperActor->GetProperty()->SetColor(1,1,0);
        renderer->AddActor(clipperActor);
        m_Controls.comboBox->insertItem(i, "Clipper " + QString::number(i+1));
        clipperActors.push_back(clipperActor);
    }//_for
    m_Controls.widget_1->GetRenderWindow()->Render();

    //Adjust controllers
    m_Controls.comboBox->setCurrentIndex(0);
    m_Controls.slider->setValue(index);
    m_Controls.slider->setEnabled(true);
    m_Controls.spinBox->setEnabled(true);
    m_Controls.comboBox->setEnabled(true);
    m_Controls.button_2->setEnabled(false);
}

void WallThicknessCalculationsClipperView::ClipperImage() {

    if (this->GetDataManagerSelection().empty()) {
        QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation to clip!");
        return;
    } else if (clipper->GetCentreLines().size() == 0) {
        QMessageBox::warning(NULL, "Attention", "Please maske sure you have computed centre lines!");
        return;
    } else if (clipper->GetCentreLinePolyPlanes().size() == 0) {
        QMessageBox::warning(NULL, "Attention", "Please maske sure you have computed clipper planes!");
        return;
    } else {

        //Check for selection of segmentation image
        QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
        mitk::DataNode::Pointer segNode = nodes.at(0);
        mitk::BaseData::Pointer data = segNode->GetData();
        if (data) {

            //Test if this data item is an image
            mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
            if (image) {

                //Check if the right segmentation
                bool morphAnalysis = false;
                if (segNode->GetName().compare(fileName.left(fileName.length()-4).toStdString()) == 0) {
                    int reply1 = QMessageBox::question(
                                NULL, "Question", "Are you sure you want to clip the same segmentation used for creating the visualised mesh?",
                                QMessageBox::Yes, QMessageBox::No);
                    if (reply1 == QMessageBox::No) {
                        QMessageBox::warning(NULL, "Attention", "Please select the correct segmentation to clip!");
                        return;
                    }//_if_no
                    int reply2 = QMessageBox::question(
                                NULL, "Question", "Do you want to save PV clipping information for morphological analysis?",
                                QMessageBox::Yes, QMessageBox::No);
                    if (reply2 == QMessageBox::Yes) {
                        QMessageBox::information(NULL, "Attention", "Labelled bloodpool image and ostia information will be saved!");
                        morphAnalysis = true;
                    }//_if_yes
                }//_if

                this->BusyCursorOn();
                this->GetDataStorage()->Remove(segNode);
                clipper->ClipVeinsImage(pickedSeedLabels, image, morphAnalysis ? true : false);
                this->BusyCursorOff();
                QMessageBox::information(NULL, "Attention", "Segmentation is now clipped!");

            } else {
                QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
                return;
            }//_image
        } else {
            QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
            return;
        }//_data
    }//_if

    //Add to storage
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(clipper->GetClippedSegImage());
    node->SetName("PVeinsCroppedImage");
    this->GetDataStorage()->Add(node);

    //Adjust controllers
    m_Controls.slider->setEnabled(false);
    m_Controls.spinBox->setEnabled(false);
    m_Controls.comboBox->setEnabled(false);
    m_Controls.button_4->setEnabled(false);
}

void WallThicknessCalculationsClipperView::CtrPlanesPlacer() {

    if (m_Controls.comboBox->count() == 0)
        return;

    //Retrieve cutting planes
    std::vector<vtkSmartPointer<vtkRegularPolygonSource>> ctrPlanes = clipper->GetCentreLinePolyPlanes();
    vtkSmartPointer<vtkRegularPolygonSource> ctrPlane = ctrPlanes.at(m_Controls.comboBox->currentIndex());

    //Parameters recalculations to adjust
    int position = m_Controls.slider->value();
    double value = m_Controls.spinBox->value();
    int indexBox = m_Controls.comboBox->currentIndex();

    clipper->SetRadiusAdjustment(value);
    clipper->CalcParamsOfPlane(ctrPlane, indexBox, position);

    ctrPlane->Update();
    m_Controls.widget_1->GetRenderWindow()->Render();
}

void WallThicknessCalculationsClipperView::CtrLinesSelector(int index) {

    if (m_Controls.comboBox->count() == 0)
        return;

    //Find the right position
    vtkSmartPointer<vtkPolyData> line = clipper->GetCentreLines().at(index)->GetOutput();
    vtkSmartPointer<vtkRegularPolygonSource> ctrPlane = clipper->GetCentreLinePolyPlanes().at(index);
    int position = line->FindPoint(ctrPlane->GetCenter());

    //Find the right radius
    vtkSmartPointer<vtkDoubleArray> radii = vtkSmartPointer<vtkDoubleArray>::New();
    radii = vtkDoubleArray::SafeDownCast(line->GetPointData()->GetArray("MaximumInscribedSphereRadius"));
    double adjust = ctrPlane->GetRadius() / radii->GetValue(position);

    //Adjust highlighted clipper
    for (unsigned int i=0; i<clipperActors.size(); i++)
        clipperActors.at(i)->GetProperty()->SetOpacity(i==(unsigned)index ? 1.0 : 0.1);

    //Set corresponding label
    int type = clipper->GetManualType().at(index);
    if (type == 0) {
        m_Controls.label->setText(" Automatic ");
        m_Controls.label->setStyleSheet("QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:yellow;}");
    } else if (type == 1) {
        m_Controls.label->setText(" Manual--1 ");
        m_Controls.label->setStyleSheet("QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:lime;}");
    } else {
        m_Controls.label->setText(" Manual--2 ");
        m_Controls.label->setStyleSheet("QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:red;}");
    }//_if_type

    //Set controllers
    m_Controls.slider->setRange(0, line->GetNumberOfPoints()-1);
    m_Controls.slider->setValue(position);
    m_Controls.spinBox->setValue(adjust);
    pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
}

void WallThicknessCalculationsClipperView::Visualiser() {

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
    surfMapper->ScalarVisibilityOff();
    vtkSmartPointer<vtkActor> surfActor = vtkSmartPointer<vtkActor>::New();
    surfActor->SetMapper(surfMapper);
    surfActor->GetProperty()->SetOpacity(1);
    renderer->AddActor(surfActor);
}

void WallThicknessCalculationsClipperView::PickCallBack() {

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

void WallThicknessCalculationsClipperView::ManualCutterCallBack() {

    if (m_Controls.label->text() != " Manual--2 ") {

        vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
        vtkSmartPointer<vtkSphereSource> glyphSource = vtkSmartPointer<vtkSphereSource>::New();
        glyph3D->SetInputData(pickedCutterSeeds);
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

        //Adjust labels
        m_Controls.label->setText(" Manual--2 ");
        m_Controls.label->setStyleSheet("QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:red;}");
    }//_if

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

    double* point = surface->GetVtkPolyData()->GetPoint(pickedSeedId);
    pickedCutterSeeds->GetPoints()->InsertNextPoint(point);
    pickedCutterSeeds->Modified();
    m_Controls.widget_1->GetRenderWindow()->Render();
}

void WallThicknessCalculationsClipperView::KeyCallBackFunc(vtkObject*, long unsigned int, void* ClientData, void*) {

    WallThicknessCalculationsClipperView* self;
    self = reinterpret_cast<WallThicknessCalculationsClipperView*>(ClientData);
    std::string key = self->interactor->GetKeySym();

    if (self->m_Controls.button_1->isEnabled()) {

        if (key == "space") {

            //Ask the labels
            self->PickCallBack();
            int dialogCode = self->inputs->exec();
            QRect screenGeometry = QApplication::desktop()->screenGeometry();
            int x = (screenGeometry.width() - self->inputs->width()) / 2;
            int y = (screenGeometry.height() - self->inputs->height()) / 2;
            self->inputs->move(x,y);

            //Act on dialog return code
            if (dialogCode == QDialog::Accepted) {

                if (self->m_Labels.radioButton_1->isChecked()) {
                    self->pickedSeedLabels.push_back(11);
                    self->m_Labels.radioButton_1->setEnabled(false);
                } else if (self->m_Labels.radioButton_2->isChecked()) {
                    self->pickedSeedLabels.push_back(12);
                    self->m_Labels.radioButton_2->setEnabled(false);
                } else if (self->m_Labels.radioButton_3->isChecked()) {
                    self->pickedSeedLabels.push_back(13);
                    self->m_Labels.radioButton_3->setEnabled(false);
                } else if (self->m_Labels.radioButton_4->isChecked()) {
                    self->pickedSeedLabels.push_back(14);
                    self->m_Labels.radioButton_4->setEnabled(false);
                } else if (self->m_Labels.radioButton_5->isChecked()) {
                    self->pickedSeedLabels.push_back(15);
                    self->m_Labels.radioButton_5->setEnabled(false);
                } else if (self->m_Labels.radioButton_6->isChecked()) {
                    self->pickedSeedLabels.push_back(16);
                    self->m_Labels.radioButton_6->setEnabled(false);
                } else if (self->m_Labels.radioButton_7->isChecked()) {
                    self->pickedSeedLabels.push_back(17);
                    self->m_Labels.radioButton_7->setEnabled(false);
                } else if (self->m_Labels.radioButton_8->isChecked()) {
                    self->pickedSeedLabels.push_back(18);
                    self->m_Labels.radioButton_8->setEnabled(false);
                } else if (self->m_Labels.radioButton_9->isChecked()) {
                    self->pickedSeedLabels.push_back(19);
                    self->m_Labels.radioButton_9->setEnabled(false);
                } else if (self->m_Labels.radioButton10->isChecked()) {
                    self->pickedSeedLabels.push_back(20);
                    self->m_Labels.radioButton10->setEnabled(false);
                } else
                    self->pickedSeedLabels.push_back(21);
                self->m_Labels.radioButton_0->setChecked(true);

            } else if (dialogCode == QDialog::Rejected) {
                self->inputs->close();
            }//_if

        } else if (key == "Delete") {

            vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
            vtkSmartPointer<vtkPoints> points = self->pickedLineSeeds->GetPoints();
            for (int i=0; i<points->GetNumberOfPoints()-1; i++)
                newPoints->InsertNextPoint(points->GetPoint(i));
            self->pickedLineSeeds->SetPoints(newPoints);
            if (self->pickedSeedLabels.empty() == false)
                self->pickedSeedLabels.pop_back();
            self->m_Controls.widget_1->GetRenderWindow()->Render();

        }//_if_space

    } else if (self->clipper->GetCentreLinePolyPlanes().size() != 0 && !self->m_Controls.button_1->isEnabled() &&
               /*self->m_Controls.button_3->isEnabled() &&*/ self->m_Controls.button_4->isEnabled()) {

        double adjustments[2];
        int idxClipper = self->m_Controls.comboBox->currentIndex();
        adjustments[0] = self->clipper->GetMClipperAngles()[idxClipper][0];
        adjustments[1] = self->clipper->GetMClipperAngles()[idxClipper][1];

        if (key == "Up") {

            self->m_Controls.label->setText(" Manual--1 ");
            QString style = "QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:lime;}";
            self->m_Controls.label->setStyleSheet(style);
            if (self->clipper->GetMClipperAngles()[idxClipper][0] >= 180.0)
                adjustments[0] = 0.0;
            else
                adjustments[0] = self->clipper->GetMClipperAngles()[idxClipper][0] + 0.1;
            self->clipperActors.at(idxClipper)->GetProperty()->SetColor(0,1,0);
            self->clipper->SetMClipperAngles(adjustments, idxClipper);
            self->CtrPlanesPlacer();

        } else if (key == "Down") {

            self->m_Controls.label->setText(" Manual--1 ");
            QString style = "QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:lime;}";
            self->m_Controls.label->setStyleSheet(style);
            if (self->clipper->GetMClipperAngles()[idxClipper][0] <= 0.0)
                adjustments[0] = 180.0;
            else
                adjustments[0] = self->clipper->GetMClipperAngles()[idxClipper][0] - 0.1;
            self->clipperActors.at(idxClipper)->GetProperty()->SetColor(0,1,0);
            self->clipper->SetMClipperAngles(adjustments, idxClipper);
            self->CtrPlanesPlacer();

        } else if (key == "Right") {

            self->m_Controls.label->setText(" Manual--1 ");
            QString style = "QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:lime;}";
            self->m_Controls.label->setStyleSheet(style);
            if (self->clipper->GetMClipperAngles()[idxClipper][1] >= 360.0)
                adjustments[1] = 0.0;
            else
                adjustments[1] = self->clipper->GetMClipperAngles()[idxClipper][1] + 0.1;
            self->clipperActors.at(idxClipper)->GetProperty()->SetColor(0,1,0);
            self->clipper->SetMClipperAngles(adjustments, idxClipper);
            self->CtrPlanesPlacer();

        } else if (key == "Left") {

            self->m_Controls.label->setText(" Manual--1 ");
            QString style = "QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:lime;}";
            self->m_Controls.label->setStyleSheet(style);
            if (self->clipper->GetMClipperAngles()[idxClipper][1] <= 0.0)
                adjustments[1] = 360.0;
            else
                adjustments[1] = self->clipper->GetMClipperAngles()[idxClipper][1] - 0.1;
            self->clipperActors.at(idxClipper)->GetProperty()->SetColor(0,1,0);
            self->clipper->SetMClipperAngles(adjustments, idxClipper);
            self->CtrPlanesPlacer();

        } else if (key == "space") {

            self->clipperActors.at(idxClipper)->GetProperty()->SetColor(1,0,0);
            self->ManualCutterCallBack();
            self->clipper->SetMClipperSeeds(self->pickedCutterSeeds, self->m_Controls.comboBox->currentIndex());
            self->CtrPlanesPlacer();

        } else if (key == "Delete") {

            vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
            vtkSmartPointer<vtkPoints> points = self->pickedCutterSeeds->GetPoints();
            for (int i=0; i<points->GetNumberOfPoints()-1; i++)
                newPoints->InsertNextPoint(points->GetPoint(i));
            self->pickedCutterSeeds->SetPoints(newPoints);
            self->clipper->SetMClipperSeeds(self->pickedCutterSeeds, self->m_Controls.comboBox->currentIndex());
            self->CtrPlanesPlacer();

        } else if (key == "a" || key == "A") {

            self->m_Controls.label->setText(" Automatic ");
            QString style = "QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:yellow;}";
            self->m_Controls.label->setStyleSheet(style);
            self->clipper->SetToAutomaticClipperMode(idxClipper);
            self->clipperActors.at(idxClipper)->GetProperty()->SetColor(1,1,0);
            self->CtrPlanesPlacer();

        } else if (key == "o" || key == "O") {

            if (self->surfActor->GetProperty()->GetOpacity() == 1.0)
                self->surfActor->GetProperty()->SetOpacity(0.5);
            else {
                self->surfActor->GetProperty()->SetOpacity(1.0);
                for (unsigned int i=0; i<self->clipperActors.size(); i++)
                    self->clipperActors.at(i)->GetProperty()->SetOpacity(1.0);
            }//_if
            self->m_Controls.widget_1->GetRenderWindow()->Render();

        }//_if_key

    } else if (!self->m_Controls.button_1->isEnabled()) {

        if (key == "r" || key =="R") {

            //Clear renderer
            self->renderer->RemoveAllViewProps();
            self->clipperActors.clear();
            self->pickedSeedLabels.clear();
            self->pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
            self->pickedSeedIds->Initialize();
            self->pickedLineSeeds = vtkSmartPointer<vtkPolyData>::New();
            self->pickedLineSeeds->Initialize();
            self->pickedLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
            self->pickedCutterSeeds = vtkSmartPointer<vtkPolyData>::New();
            self->pickedCutterSeeds->Initialize();
            self->pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
            self->Visualiser();

            //Reset clipper
            if (self->m_Controls.button_2->isEnabled() && self->m_Controls.button_4->isEnabled())
                self->clipper = std::unique_ptr<CemrgAtriaClipper>(new CemrgAtriaClipper(self->directory, self->surface));
            else
                self->clipper->ResetCtrLinesClippingPlanes();

            //Reset controls
            self->m_Controls.comboBox->clear();
            self->m_Controls.slider->setRange(0,2);
            self->m_Controls.slider->setValue(0);
            self->m_Controls.spinBox->setValue(2.0);
            self->m_Controls.slider->setEnabled(false);
            self->m_Controls.spinBox->setEnabled(false);
            self->m_Controls.comboBox->setEnabled(false);
            self->m_Controls.button_1->setEnabled(true);
            self->m_Controls.button_2->setEnabled(true);
            self->m_Controls.button_4->setEnabled(true);
            self->m_Controls.checkBox->setEnabled(true);
            self->m_Labels.radioButton_1->setEnabled(true);
            self->m_Labels.radioButton_2->setEnabled(true);
            self->m_Labels.radioButton_3->setEnabled(true);
            self->m_Labels.radioButton_4->setEnabled(true);
            self->m_Labels.radioButton_5->setEnabled(true);
            self->m_Labels.radioButton_6->setEnabled(true);
            self->m_Labels.radioButton_7->setEnabled(true);
            self->m_Labels.radioButton_8->setEnabled(true);
            self->m_Labels.radioButton_9->setEnabled(true);
            self->m_Labels.radioButton10->setEnabled(true);

            //Setup shortcuts
            vtkSmartPointer<vtkTextActor> txtActor = vtkSmartPointer<vtkTextActor>::New();
            std::string shortcuts = "R: reset centrelines\nSpace: add seed point\nDelete: remove seed point";
            txtActor->SetInput(shortcuts.c_str());
            txtActor->GetTextProperty()->SetFontSize(14);
            txtActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
            self->renderer->AddActor2D(txtActor);

        }//_if
    }//_if_main
}
