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
 * Motion Quantification
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
#include <berryFileEditorInput.h>

// Qmitk
#include <QmitkPlotWidget.h>
#include <QmitkIOUtil.h>
#include <mitkCoreObjectFactory.h>
#include <mitkDataNode.h>
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkPlanarCircle.h>
#include "MmcwViewPlot.h"

// VTK
#include <vtkLineSource.h>
#include <vtkSectorSource.h>
#include <vtkRegularPolygonSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkColorTransferFunction.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkScalarBarActor.h>
#include <vtkCellData.h>

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>
#include <QSignalMapper>

QString MmcwViewPlot::directory;
int MmcwViewPlot::noFrames = 10;
int MmcwViewPlot::smoothness = 1;
const std::string MmcwViewPlot::VIEW_ID = "org.mitk.views.mmcwplot";

void MmcwViewPlot::SetFocus() {
}

void MmcwViewPlot::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, &QPushButton::clicked, this, &MmcwViewPlot::PlotData);
    connect(m_Controls.button_2, &QPushButton::clicked, this, &MmcwViewPlot::FilePlot);
    connect(m_Controls.button_3, &QPushButton::clicked, this, &MmcwViewPlot::BullPlot);
    //connect(m_Controls.horizontalSlider, SIGNAL(valueChanged(int)), this, SLOT(ColourAHASegments(int)));
    connect(m_Controls.horizontalSlider, &QSlider::valueChanged, this, &MmcwViewPlot::ColourAHASegments);

    //Adjust controllers
    m_Controls.lineEdit_F->setPlaceholderText("No Frames (default = " + QString::number(noFrames) + ")");
    m_Controls.comboBox_S->setCurrentIndex(smoothness < 5 ? smoothness-1 : 2);
    m_Controls.horizontalSlider->setMaximum(noFrames*smoothness);
    cardiCycle = 0;

    //AHA bullseye plot
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow =
            vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    m_Controls.widget_1->SetRenderWindow(renderWindow);

    AHA_renderer = vtkSmartPointer<vtkRenderer>::New();
    AHA_renderer->SetBackground(0,0,0);
    m_Controls.horizontalSlider->setMaximum(noFrames*smoothness);
    m_Controls.widget_1->GetRenderWindow()->AddRenderer(AHA_renderer);
    AHA_interactor = m_Controls.widget_1->GetRenderWindow()->GetInteractor();
    AHA_interactor->RemoveObservers(vtkCommand::LeftButtonPressEvent);
    AHA_interactor->RemoveObservers(vtkCommand::LeftButtonReleaseEvent);
    AHA_interactor->RemoveObservers(vtkCommand::RightButtonPressEvent);
    AHA_interactor->RemoveObservers(vtkCommand::RightButtonReleaseEvent);
    //AHA_interactor->RemoveObservers(vtkCommand::MouseWheelForwardEvent);
    //AHA_interactor->RemoveObservers(vtkCommand::MouseWheelBackwardEvent);
    AHA_interactor->Start();

    //Plot section
    m_Controls.widget_2->SetPlotTitle("AHA Curves Plot");
    m_Controls.widget_2->SetAxisTitle(QwtPlot::xBottom, "Time");
    m_Controls.widget_2->SetAxisTitle(QwtPlot::yLeft, "Value");

    //Setup camera
    AHA_camera = vtkSmartPointer<vtkCamera>::New();
    AHA_renderer->SetActiveCamera(AHA_camera);

    //Setup AHA lookup table
    AHA[1]=13; AHA[5]=10; AHA[9 ]=8; AHA[13]=6;
    AHA[2]=14; AHA[6]=11; AHA[10]=9; AHA[14]=1;
    AHA[3]=15; AHA[7]=12; AHA[11]=4; AHA[15]=2;
    AHA[4]=16; AHA[8]= 7; AHA[12]=5; AHA[16]=3;
}

void MmcwViewPlot::SetDirectory(const QString directory) {

    MmcwViewPlot::directory = directory;
}

void MmcwViewPlot::SetNoFrames(int frames, int smoothness) {

    MmcwViewPlot::noFrames = frames;
    MmcwViewPlot::smoothness = smoothness;
}

void MmcwViewPlot::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void MmcwViewPlot::PlotData() {

    //Check for selection of landmarks
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(NULL, "Attention", "Please select landmarks from the Data Manager to continue!");
        return;
    }

    //Ask the user for a dir to locate data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }

    //Frames and smoothness adjustments
    bool ok1;
    int frames = m_Controls.lineEdit_F->text().toInt(&ok1);
    smoothness = m_Controls.comboBox_S->currentIndex() < 2 ? m_Controls.comboBox_S->currentIndex()+1 : 5;
    noFrames = (ok1) ? frames : noFrames;
    m_Controls.horizontalSlider->setMaximum(noFrames*smoothness);

    //Find the reference mesh
    bool ok2;
    int refMshNo = QInputDialog::getInt(NULL, tr("Reference Mesh"), tr("Number:"), 0, 0, noFrames*smoothness, 1, &ok2);
    if (!ok2) {
        QMessageBox::warning(NULL, "Attention", "Have you completed the tracking step?");
        return;
    }

    //Segments sections ratio
    int bas = m_Controls.lineEdit_1->text().toInt();
    int mid = m_Controls.lineEdit_2->text().toInt();
    int api = m_Controls.lineEdit_3->text().toInt();
    if (std::abs(bas+mid+api - 100) > 1.0) {
        QMessageBox::warning(NULL, "Attention", "Revert to a default ratio for basal, mid, and apical segments!");
        bas = 33; mid = 33; api = 33;
    }

    this->BusyCursorOn();
    AHA_renderer->RemoveAllViewProps();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(1);

    //Define the reference mesh
    mitk::Surface::Pointer refSurf;
    int segRatios[3] = {bas, mid, api};
    int pacingSegRatios[3] = {17, 33, 55};
    // This give the AHA region from 50 - 50+1/3 in middle - Ie putting pacing site at 2/3 height
    strain = std::unique_ptr<CemrgStrains>(new CemrgStrains(directory, refMshNo));
    mitk::DataNode::Pointer lmNode = nodes.front();
    if (!dynamic_cast<mitk::PointSet*>(lmNode->GetData())) {
        QMessageBox::warning(NULL, "Attention", "Please select landmarks from the Data Manager to continue!");
        this->BusyCursorOff();
        mitk::ProgressBar::GetInstance()->Progress();
        return;
    }

    //Calculate y values of the plots
    flatPlotScalars.clear();
    plotValueVectors.clear();
    std::string plotType = m_Controls.comboBox->currentText().toStdString();

    if (plotType.compare("Area Change") == 0) {
        refSurf = strain->ReferenceAHA(lmNode, segRatios, false);
        for (int i=0; i<noFrames*smoothness; i++) {
            plotValueVectors.push_back(strain->CalculateSqzPlot(i));
            vtkSmartPointer<vtkFloatArray> holder = vtkSmartPointer<vtkFloatArray>::New();
            holder->DeepCopy(strain->GetFlatSurfScalars());
            flatPlotScalars.push_back(holder);
        }
    } else if (plotType.compare("Circumferential Small Strain") == 0) {
        refSurf = strain->ReferenceAHA(lmNode, segRatios, false);
        for (int i=0; i<noFrames*smoothness; i++) {
            plotValueVectors.push_back(strain->CalculateStrainsPlot(i, lmNode, 1));
            vtkSmartPointer<vtkFloatArray> holder = vtkSmartPointer<vtkFloatArray>::New();
            holder->DeepCopy(strain->GetFlatSurfScalars());
            flatPlotScalars.push_back(holder);
        }
    } else if (plotType.compare("Circumferential Large Strain") == 0) {
        refSurf = strain->ReferenceAHA(lmNode, segRatios, false);
        for (int i=0; i<noFrames*smoothness; i++) {
            plotValueVectors.push_back(strain->CalculateStrainsPlot(i, lmNode, 3));
            vtkSmartPointer<vtkFloatArray> holder = vtkSmartPointer<vtkFloatArray>::New();
            holder->DeepCopy(strain->GetFlatSurfScalars());
            flatPlotScalars.push_back(holder);
        }
    } else if (plotType.compare("Longitudinal Small Strain") == 0) {
        refSurf = strain->ReferenceAHA(lmNode, segRatios, false);
        for (int i=0; i<noFrames*smoothness; i++) {
            plotValueVectors.push_back(strain->CalculateStrainsPlot(i, lmNode, 2));
            vtkSmartPointer<vtkFloatArray> holder = vtkSmartPointer<vtkFloatArray>::New();
            holder->DeepCopy(strain->GetFlatSurfScalars());
            flatPlotScalars.push_back(holder);
        }
    } else if (plotType.compare("Longitudinal Large Strain") == 0) {
        refSurf = strain->ReferenceAHA(lmNode, segRatios, false);
        for (int i=0; i<noFrames*smoothness; i++) {
            plotValueVectors.push_back(strain->CalculateStrainsPlot(i, lmNode, 4));
            vtkSmartPointer<vtkFloatArray> holder = vtkSmartPointer<vtkFloatArray>::New();
            holder->DeepCopy(strain->GetFlatSurfScalars());
            flatPlotScalars.push_back(holder);
        }
    } else if (plotType.compare("Pacing site Squeez") == 0) {
        refSurf = strain->ReferenceAHA(lmNode, pacingSegRatios, true);
        for (int i=0; i<noFrames*smoothness; i++) {
            plotValueVectors.push_back(strain->CalculateSqzPlot(i));
            vtkSmartPointer<vtkFloatArray> holder = vtkSmartPointer<vtkFloatArray>::New();
            holder->DeepCopy(strain->GetFlatSurfScalars());
            flatPlotScalars.push_back(holder);
        }
    }//_if

    //Visualise AHA plots
    HandleBullPlot(false);
    ColourAHASegments(m_Controls.horizontalSlider->value());

    //Visualise the curves plot
    HandleCurvPlot();

    //Visualise the reference mesh
    mitk::DataStorage::SetOfObjects::ConstPointer sob = this->GetDataStorage()->GetAll();
    for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt) {
        if (nodeIt->Value()->GetName().compare("Segmented Mesh") == 0)
            this->GetDataStorage()->Remove(nodeIt->Value());
        if (nodeIt->Value()->GetName().compare("Reference Mesh") == 0)
            this->GetDataStorage()->Remove(nodeIt->Value());
        if (nodeIt->Value()->GetName().compare("Guideline 1") == 0)
            this->GetDataStorage()->Remove(nodeIt->Value());
        if (nodeIt->Value()->GetName().compare("Guideline 2") == 0)
            this->GetDataStorage()->Remove(nodeIt->Value());
        if (nodeIt->Value()->GetName().compare("Guideline 3") == 0)
            this->GetDataStorage()->Remove(nodeIt->Value());
    }//_for
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetProperty("scalar visibility", mitk::BoolProperty::New(true));
    node->SetName("Reference Mesh");
    node->SetData(refSurf);
    this->GetDataStorage()->Add(node);

    //Visualise reference mesh guidelines
    std::vector<mitk::Surface::Pointer> guidelines = strain->ReferenceGuideLines(lmNode);
    for (int i=0; i<3; i++) {
        mitk::DataNode::Pointer gNode = mitk::DataNode::New();
        gNode->SetName("Guideline "+ std::to_string(i+1));
        gNode->SetData(guidelines.at(i));
        this->GetDataStorage()->Add(gNode, node);
    }

    this->BusyCursorOff();
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
    mitk::ProgressBar::GetInstance()->Progress();
}

void MmcwViewPlot::BullPlot() {

    if (!strain) {
        //if plot values have not been calculated
        return;
    } else if (m_Controls.button_3->isChecked() == false) {
        HandleBullPlot(false);
        m_Controls.button_3->setText("Global");
    } else {
        HandleBullPlot(true);
        m_Controls.button_3->setText("Local");
    }//_if
    ColourAHASegments(m_Controls.horizontalSlider->value());
}

void MmcwViewPlot::FilePlot() {

    //Check anything to plot
    if (plotValueVectors.size() == 0) {
        QMessageBox::warning(NULL, "Attention", "No plot to save!");
        return;
    }

    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }

    //Write to CSV or VTK file
    if (m_Controls.button_3->isChecked() == false)
        WritePlotToCSV(directory);
    else
        WritePlotToVTK(directory);
}

void MmcwViewPlot::ColourAHASegments(int /*value*/) {

    if (plotValueVectors.size() != 0) {
        if (m_Controls.button_3->isChecked() == false)
            HandleBullPlot(false);
        else
            HandleBullPlot(true);

        //SDI information
        DrawAHATextInfo();

        //Render the window
        AHA_camera->SetPosition(0,0,m_Controls.button_3->isChecked()?300:12);
        AHA_renderer->GetRenderWindow()->Render();
    }//_if
}

/**************************************************************************************************
 *************** PRIVATE FUNCTIONS ****************************************************************
 **************************************************************************************************/

void MmcwViewPlot::HandleBullPlot(bool global) {

    //Clean up the renderer
    AHA_renderer->RemoveAllViewProps();

    if (global == false) {

        //Setup the surface
        std::vector<double> ranges;
        mitk::Surface::Pointer surface = strain->FlattenedAHA();
        for (int i=0; i<noFrames*smoothness; i++) {
            surface->GetVtkPolyData()->GetCellData()->SetScalars(flatPlotScalars.at(i));
            ranges.push_back(surface->GetVtkPolyData()->GetScalarRange()[0]);
            ranges.push_back(surface->GetVtkPolyData()->GetScalarRange()[1]);
        }

        //Setup lookup table
        double range[2];
        range[0] = *std::min_element(ranges.begin(), ranges.end());
        range[1] = *std::max_element(ranges.begin(), ranges.end());
        vtkSmartPointer<vtkColorTransferFunction> lut = GetLookupTable(range);

        //Setup scalar bar
        vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
        scalarBar->SetLookupTable(lut);
        scalarBar->SetNumberOfLabels(2);
        scalarBar->SetPosition(0.9,0.1);
        scalarBar->SetWidth(0.1);
        scalarBar->SetTitle(" ");
        AHA_renderer->AddActor2D(scalarBar);

        //Setup AHA segments
        int frame = (m_Controls.horizontalSlider->value() == noFrames*smoothness) ? 0 : m_Controls.horizontalSlider->value();
        DrawAHASegments(frame, range);
        DrawAHALines();

    } else {

        //Setup the surface
        std::vector<double> ranges;
        mitk::Surface::Pointer surface = strain->FlattenedAHA();
        for (int i=0; i<noFrames*smoothness; i++) {
            surface->GetVtkPolyData()->GetCellData()->SetScalars(flatPlotScalars.at(i));
            ranges.push_back(surface->GetVtkPolyData()->GetScalarRange()[0]);
            ranges.push_back(surface->GetVtkPolyData()->GetScalarRange()[1]);
        }
        int frame = (m_Controls.horizontalSlider->value() == noFrames*smoothness) ? 0 : m_Controls.horizontalSlider->value();
        surface->GetVtkPolyData()->GetCellData()->SetScalars(flatPlotScalars.at(frame));

        //Setup lookup table
        double range[2];
        range[0] = *std::min_element(ranges.begin(), ranges.end());
        range[1] = *std::max_element(ranges.begin(), ranges.end());
        vtkSmartPointer<vtkColorTransferFunction> lut = GetLookupTable(range);

        //Setup scalar bar
        vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
        scalarBar->SetLookupTable(lut);
        scalarBar->SetNumberOfLabels(2);
        scalarBar->SetPosition(0.9,0.1);
        scalarBar->SetWidth(0.1);
        scalarBar->SetTitle(" ");
        AHA_renderer->AddActor2D(scalarBar);

        //Create a mapper and actor
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(surface->GetVtkPolyData());
        mapper->ScalarVisibilityOn();
        mapper->SetScalarModeToUseCellData();
        mapper->SetScalarRange(range[0], range[1]);
        mapper->SetLookupTable(lut);
        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->LightingOff();
        AHA_renderer->AddActor2D(actor);
    }//_if
}

void MmcwViewPlot::HandleCurvPlot() {

    int curveId = 0;
    bool nonComputable = false;
    m_Controls.widget_2->Clear();
    QmitkPlotWidget::DataVector xValues(noFrames*smoothness+1,0);
    QmitkPlotWidget::DataVector yValues(noFrames*smoothness+1,0);

    for (int i=0; i<16; i++) {
        //Order x and y values
        for (int j=0; j<noFrames*smoothness+1; j++) {
            xValues[j] = j;
            if (j<noFrames*smoothness)
                yValues[j] = plotValueVectors[j][i];
            else
                yValues[j] = plotValueVectors[0][i];
            if (std::isnan(yValues[j]))
                nonComputable = true;
        }//_for

        //Setup the curve
        int label = i + 1;
        std::string title = "Segment " + QString::number(label).toStdString();
        std::vector<float> colour = strain->GetAHAColour(label);
        QColor qColour(colour[0], colour[1], colour[2]);
        legend = std::unique_ptr<QwtLegend>(new QwtLegend());
        curveId = m_Controls.widget_2->InsertCurve(std::to_string(label).c_str());

        m_Controls.widget_2->SetCurveData(curveId, xValues, yValues);
        m_Controls.widget_2->SetCurvePen(curveId, QPen(qColour));
        m_Controls.widget_2->SetCurveTitle(curveId, title.c_str());
        m_Controls.widget_2->SetPlotTitle("AHA Curves Plot");
        m_Controls.widget_2->SetLegend(legend.get(), QwtPlot::RightLegend, 0.5);
        m_Controls.widget_2->Replot();
    }//_for

    QString msg =
            "Curves were not computed. Check the order of selected landmarks! "
            "For instance, order of two points on RV cusps.";
    if (nonComputable)
        QMessageBox::warning(NULL, "Attention", msg);
}

void MmcwViewPlot::DrawAHALines() {

    //Draw lines
    const double z = 0.0001;
    //Create two points, P0 and P1
    std::vector<double> p0, p1;

    for (int i=0; i<8; i++) {

        vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
        if (i==0) {
            p0 = {-.70,-.70,z}; p1 = {0.70,0.70,z};
        } else if (i==1) {
            p0 = {-.70,0.70,z}; p1 = {0.70,-.70,z};
        } else if (i==2) {
            p0 = {1.00,0.00,z}; p1 = {3.00,0.00,z};
        } else if (i==3) {
            p0 = {-1.0,0.00,z}; p1 = {-3.0,0.00,z};
        } else if (i==4) {
            p0 = {0.50,0.86,z}; p1 = {1.50,2.59,z};
        } else if (i==5) {
            p0 = {0.50,-.86,z}; p1 = {1.5,-2.59,z};
        } else if (i==6) {
            p0 = {-.50,0.86,z}; p1 = {-1.5,2.59,z};
        } else {
            p0 = {-.5,-.86,z}; p1 = {-1.5,-2.59,z};
        }//_if

        lineSource->SetPoint1(p0[0], p0[1], p0[2]);
        lineSource->SetPoint2(p1[0], p1[1], p1[2]);
        lineSource->Update();

        //Create a mapper and actor
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(lineSource->GetOutputPort());
        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);

        //Setup the colour
        actor->GetProperty()->SetLineWidth(3);
        actor->GetProperty()->SetColor(0.0,0.0,0.0);

        //Add to renderer
        AHA_renderer->AddActor2D(actor);
    }

    //Draw Circles
    for (int i=0; i<2; i++) {

        vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
        polygonSource->GeneratePolygonOff();
        polygonSource->SetNumberOfSides(50);
        polygonSource->SetRadius(i+1);
        polygonSource->SetCenter(0.0,0.0,0.0);

        //Create a mapper and actor
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(polygonSource->GetOutputPort());
        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);

        //Setup the colour
        actor->GetProperty()->SetLineWidth(3);
        actor->GetProperty()->SetColor(0.0,0.0,0.0);

        //Add to renderer
        AHA_renderer->AddActor2D(actor);
    }//_for
}

void MmcwViewPlot::DrawAHASegments(int frame, double* range) {

    //Setup lookup table
    double angle, radii, radio, inf = 0.0001;
    vtkSmartPointer<vtkColorTransferFunction> lut = GetLookupTable(range);

    //Create Segments Layers
    for (int i=0; i<16; i++) {

        vtkSmartPointer<vtkSectorSource> segmentSource = vtkSmartPointer<vtkSectorSource>::New();
        if (i<4) {
            angle = 90; radii = inf; radio = 1.0;
        } else if (i<10) {
            angle = 60; radii = 1.0; radio = 2.0;
        } else {
            angle = 60; radii = 2.0; radio = 3.0;
        }//_if
        segmentSource->SetStartAngle(i*angle);
        segmentSource->SetEndAngle((i+1)*angle);
        segmentSource->SetInnerRadius(radii);
        segmentSource->SetOuterRadius(radio);

        //Create a mapper and actor
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(segmentSource->GetOutputPort());
        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);

        //Adjust rotation
        if (i<4)
            actor->RotateZ(45.0);

        //Setup colours
        double colour[3];
        double value = plotValueVectors[frame][AHA.find(i+1)->second-1];
        lut->GetColor(value, colour);
        actor->GetProperty()->SetColor(colour[0], colour[1], colour[2]);

        //Add labels to actors
        float offst = (AHA.find(i+1)->second > 9) ? .15 : .1;
        double* pos = vtkProp3D::SafeDownCast(actor)->GetCenter();
        vtkSmartPointer<vtkTextActor> textActor = vtkSmartPointer<vtkTextActor>::New();
        textActor->GetPositionCoordinate()->SetCoordinateSystemToWorld();
        textActor->SetInput(std::to_string(AHA.find(i+1)->second).c_str());
        textActor->SetPosition(pos[0]-offst, pos[1]-offst);
        textActor->GetTextProperty()->SetFontSize(14);
        textActor->GetTextProperty()->SetColor(0.0, 0.0, 0.0);

        //Add to renderer
        AHA_renderer->AddActor2D(actor);
        AHA_renderer->AddActor2D(textActor);
    }//_for
}

void MmcwViewPlot::DrawAHATextInfo() {

    bool ok = true;
    if (cardiCycle == 0)
        cardiCycle = QInputDialog::getInt(NULL, tr("Cycle Length in ms"), tr("Value:"), 1000, 1, 2000, 1, &ok);
    if (ok) {
        double SDI = strain->CalculateSDI(plotValueVectors, cardiCycle, noFrames*smoothness);
        std::ostringstream os;
        os << std::fixed << std::setprecision(2) << SDI;
        std::string output = "SDI: " + os.str() + "%";
        //Setup the text and add it to the renderer
        vtkSmartPointer<vtkTextActor> textActor = vtkSmartPointer<vtkTextActor>::New();
        textActor->SetInput(output.c_str());
        textActor->GetTextProperty()->SetFontSize(20);
        textActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
        AHA_renderer->AddActor2D(textActor);
    } else
        cardiCycle = 0;
}

void MmcwViewPlot::WritePlotToCSV(QString dir) {

    if (plotValueVectors.size() != 0) {

        bool ok;
        QString fileName = "Plot.csv";
        if (m_Controls.comboBox->currentText().startsWith("A"))
            fileName = "SQZ.csv";
        else if (m_Controls.comboBox->currentText().startsWith("C"))
            fileName = "CRC.csv";
        else if (m_Controls.comboBox->currentText().startsWith("L"))
            fileName = "LNG.csv";
        else if (m_Controls.comboBox->currentText().startsWith("P"))
            fileName = "SQZ_at_pacing_0.66.csv";
        fileName = QInputDialog::getText(NULL, tr("Save As"), tr("File Name:"), QLineEdit::Normal, fileName, &ok);
        if (ok && !fileName.isEmpty() && fileName.endsWith(".csv")) {

            ofstream file;
            file.open(dir.toStdString() + mitk::IOUtil::GetDirectorySeparator() + fileName.toStdString());
            std::vector<double> values;
            for (int i=0; i<16; i++) {
                for (int j=0; j<noFrames*smoothness; j++)
                    values.push_back(plotValueVectors[j][i]);
                //Append the curve to the file
                for (size_t z=0; z<values.size(); z++) {
                    file << values.at(z);
                    if (z == values.size()-1) file << endl;
                    else file << ",";
                }
                values.clear();
            }//_for
            file.close();
            QMessageBox::information(NULL, "Attention", "The plot was saved to a CSV file.");

        } else {
            QMessageBox::warning(NULL, "Attention", "Please type a file name with the right extension (i.e. .csv)!");
            return;
        }//_fileName
    }//_if
}

void MmcwViewPlot::WritePlotToVTK(QString dir) {

    mitk::Surface::Pointer surface = strain->FlattenedAHA();
    int frame = (m_Controls.horizontalSlider->value() == noFrames*smoothness) ? 0 : m_Controls.horizontalSlider->value();
    surface->GetVtkPolyData()->GetCellData()->SetScalars(flatPlotScalars.at(frame));
    QString fileName = dir + mitk::IOUtil::GetDirectorySeparator() + QString::number(frame) + ".vtk";
    mitk::IOUtil::Save(surface, fileName.toStdString());
    QMessageBox::information(NULL, "Attention", "The plot was saved to a VTK file.");
}

vtkSmartPointer<vtkColorTransferFunction> MmcwViewPlot::GetLookupTable(double* range) {

    double middlePt = 0.0; //(range[0] + range[1]) / 2.0;
    vtkSmartPointer<vtkColorTransferFunction> lut = vtkSmartPointer<vtkColorTransferFunction>::New();
    lut->SetColorSpaceToRGB();
    lut->AddRGBPoint(range[0], 0,0,1); //blue
    lut->AddRGBPoint(middlePt, 1,1,1); //white
    lut->AddRGBPoint(range[1], 1,0,0); //red
    lut->SetScaleToLinear();
    return lut;
}
