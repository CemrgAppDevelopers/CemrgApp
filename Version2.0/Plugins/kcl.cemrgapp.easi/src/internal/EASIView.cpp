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
 * Eikonal Activation Simulation (EASI) Plugin for MITK
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
#include <berryFileEditorInput.h>

// Qmitk
#include <QmitkIOUtil.h>
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkBoundingObject.h>
#include <mitkCuboid.h>
#include <mitkEllipsoid.h>
#include <mitkAffineImageCropperInteractor.h>
#include <mitkIDataStorageService.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkDataStorageEditorInput.h>
#include <mitkUnstructuredGrid.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImageCast.h>
#include <mitkScaleOperation.h>
#include <mitkInteractionConst.h>
#include <mitkImage.h>
#include "kcl_cemrgapp_easi_Activator.h"
#include "EASIView.h"
#include "FibresView.h"

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>

// VTK
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

// CemrgAppModule
#include <CemrgCommandLine.h>
#include <CemrgCommonUtils.h>
// #include "CemrgTests.cpp" // Might have to be accessed directly to a path in the developer's specific workstation

// Generic
#include <usModuleRegistry.h>
#include <numeric>
#include <iostream>

const std::string EASIView::VIEW_ID = "org.mitk.views.easi";

void EASIView::SetFocus() {

    m_Controls.button_1->setFocus();
}

void EASIView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));
    connect(m_Controls.button_3, SIGNAL(clicked()), this, SLOT(SegmentIMGS()));
    connect(m_Controls.button_4, SIGNAL(clicked()), this, SLOT(CreateMesh()));
    connect(m_Controls.btn_4_bis, SIGNAL(clicked()), this, SLOT(LoadFibresView()));
    connect(m_Controls.button_5, SIGNAL(clicked()), this, SLOT(ActivationSites()));
    connect(m_Controls.button_6, SIGNAL(clicked()), this, SLOT(Simulation()));
    connect(m_Controls.button_l, SIGNAL(clicked()), this, SLOT(LoadMesh()));
    connect(m_Controls.button_r, SIGNAL(clicked()), this, SLOT(Reset()));

    //Set visibility of buttons
    m_Controls.button_2_1->setVisible(false);
    m_Controls.button_2_2->setVisible(false);
    m_Controls.button_2_3->setVisible(false);
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));
    connect(m_Controls.button_2_2, SIGNAL(clicked()), this, SLOT(CropinIMGS()));
    connect(m_Controls.button_2_3, SIGNAL(clicked()), this, SLOT(ResampIMGS()));
    m_Controls.button_5_1->setVisible(false);
    connect(m_Controls.button_5_1, SIGNAL(clicked()), this, SLOT(ConfrmSITE()));
}

void EASIView::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void EASIView::LoadDICOM() {

    //Use MITK DICOM editor
    QString editor_id = "org.mitk.editors.dicomeditor";
    berry::IEditorInput::Pointer input(new berry::FileEditorInput(QString()));
    this->GetSite()->GetPage()->OpenEditor(input, editor_id);
}

void EASIView::ProcessIMGS() {

    //Toggle visibility of buttons
    if (m_Controls.button_2_1->isVisible()) {
        m_Controls.button_2_1->setVisible(false);
        m_Controls.button_2_2->setVisible(false);
        m_Controls.button_2_3->setVisible(false);
    } else {
        m_Controls.button_2_1->setVisible(true);
        m_Controls.button_2_2->setVisible(true);
        m_Controls.button_2_3->setVisible(true);
    }
}

void EASIView::ConvertNII() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 10) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please load and select 10 images from the Data Manager before starting this step!");
        return;
    }//_if

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

    //Order dicoms based on their cycle stages
    std::vector<int> indexNodes;
    std::string seriesDescription;
    foreach (mitk::DataNode::Pointer node, nodes) {
        node->GetData()->GetPropertyList()->GetStringProperty("dicom.series.SeriesDescription", seriesDescription);
        if (seriesDescription.find("90.0%")      != seriesDescription.npos) indexNodes.push_back(9);
        else if (seriesDescription.find("80.0%") != seriesDescription.npos) indexNodes.push_back(8);
        else if (seriesDescription.find("70.0%") != seriesDescription.npos) indexNodes.push_back(7);
        else if (seriesDescription.find("60.0%") != seriesDescription.npos) indexNodes.push_back(6);
        else if (seriesDescription.find("50.0%") != seriesDescription.npos) indexNodes.push_back(5);
        else if (seriesDescription.find("40.0%") != seriesDescription.npos) indexNodes.push_back(4);
        else if (seriesDescription.find("30.0%") != seriesDescription.npos) indexNodes.push_back(3);
        else if (seriesDescription.find("20.0%") != seriesDescription.npos) indexNodes.push_back(2);
        else if (seriesDescription.find("10.0%") != seriesDescription.npos) indexNodes.push_back(1);
        else if (seriesDescription.find("0.0%")  != seriesDescription.npos) indexNodes.push_back(0);
    }//_for

    //Sort indexes based on comparing values
    std::vector<int> index(indexNodes.size());
    std::iota(index.begin(), index.end(), 0);
    std::sort(index.begin(), index.end(), [&](int i1, int i2) {return indexNodes[i1]<indexNodes[i2];});
    //Warning for cases when order is not found
    size_t length1 = nodes.size();
    size_t length2 = indexNodes.size();
    if (length1 != length2) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Cannot find the order of images automatically. Revert to user order and selections in the data manager!");
        index.resize(nodes.size());
        std::iota(index.begin(), index.end(), 0);
    }//_if

    //Convert to Nifti
    int ctr = 0;
    QString path;
    bool successfulNitfi;

    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(index.size());
    foreach (int idx, index) {
        path = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + QString::number(ctr++) + ".nii";
        successfulNitfi = CemrgCommonUtils::ConvertToNifti(nodes.at(idx)->GetData(), path);
        if (successfulNitfi) {
            this->GetDataStorage()->Remove(nodes.at(idx));
        } else {
            mitk::ProgressBar::GetInstance()->Progress(index.size());
            return;
        }//_if
        mitk::ProgressBar::GetInstance()->Progress();
    }//for
    nodes.clear();
    this->BusyCursorOff();

    //Load first item
    ctr = 0;
    path = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + QString::number(ctr) + ".nii";
    mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
}

void EASIView::CropinIMGS() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please select an image from the Data Manager to perform cropping!");
        return;
    }//_if

    //Check to cut now or not
    if (m_Controls.button_2_2->text() == QString::fromStdString("Are you done?")) {

        QString path;
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

        //Cut selected image
        this->BusyCursorOn();
        mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
        mitk::Image::Pointer outputImage = CemrgCommonUtils::CropImage();
        path = directory + mitk::IOUtil::GetDirectorySeparator() + CemrgCommonUtils::GetImageNode()->GetName().c_str() + ".nii";
        mitk::IOUtil::Save(outputImage, path.toStdString());
        mitk::ProgressBar::GetInstance()->Progress();
        this->BusyCursorOff();

        //Update datastorage
        CemrgCommonUtils::AddToStorage(outputImage, CemrgCommonUtils::GetImageNode()->GetName(), this->GetDataStorage());
        this->GetDataStorage()->Remove(CemrgCommonUtils::GetImageNode());
        this->GetDataStorage()->Remove(CemrgCommonUtils::GetCuttingNode());
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

        m_Controls.button_2_2->setText("Crop Images");
        return;
    }//_if

    //Prepare cutting cuboid
    mitk::Cuboid::Pointer cuttingObject = mitk::Cuboid::New();
    mitk::DataNode::Pointer cuttingNode = mitk::DataNode::New();
    cuttingNode->SetData(cuttingObject);
    cuttingNode->SetProperty("opacity", mitk::FloatProperty::New(0.4));
    cuttingNode->SetProperty("color", mitk::ColorProperty::New(1.0, 1.0, 0.0));
    cuttingNode->SetProperty("name", mitk::StringProperty::New("Cropper"));
    this->GetDataStorage()->Add(cuttingNode);

    //Mouse interactions
    mitk::AffineImageCropperInteractor::Pointer affineDataInteractor = mitk::AffineImageCropperInteractor::New();
    affineDataInteractor->LoadStateMachine("ClippingPlaneInteraction3D.xml", us::ModuleRegistry::GetModule("MitkDataTypesExt"));
    affineDataInteractor->SetEventConfig("CropperDeformationConfig.xml", us::ModuleRegistry::GetModule("MitkDataTypesExt"));
    affineDataInteractor->SetDataNode(cuttingNode);
    cuttingNode->SetBoolProperty("pickable", true);

    //Fit the cuboid to the image
    mitk::Image::Pointer imageToCut;
    mitk::BoundingObject::Pointer cuttingCube;
    mitk::DataNode::Pointer imageNode = nodes.at(0);
    mitk::BaseData::Pointer data = imageNode->GetData();
    cuttingCube = dynamic_cast<mitk::BoundingObject*>(cuttingNode->GetData());
    if (data) {
        //Test if this data item is an image
        imageToCut = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (imageToCut)
            cuttingCube->FitGeometry(imageToCut->GetGeometry());
        else return;
    } else return;

    //To be used for actual cutting
    CemrgCommonUtils::SetImageToCut(imageToCut);
    CemrgCommonUtils::SetCuttingCube(cuttingCube);
    CemrgCommonUtils::SetImageNode(imageNode);
    CemrgCommonUtils::SetCuttingNode(cuttingNode);
    m_Controls.button_2_2->setText("Are you done?");
}

void EASIView::ResampIMGS() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please select an image from the Data Manager to perform downsampling!");
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

    //Find the selected node
    QString path;
    mitk::DataNode::Pointer imgNode = nodes.at(0);
    mitk::BaseData::Pointer data = imgNode->GetData();
    if (data) {

        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {

            bool ok;
            int factor = QInputDialog::getInt(NULL, tr("Downsampling"), tr("By factor of:"), 3, 1, 5, 1, &ok);
            if (ok) {

                //Downsample selected image
                this->BusyCursorOn();
                mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
                mitk::Image::Pointer outputImage = CemrgCommonUtils::Downsample(image, factor);
                path = directory + mitk::IOUtil::GetDirectorySeparator() + imgNode->GetName().c_str() + ".nii";
                mitk::IOUtil::Save(outputImage, path.toStdString());
                mitk::ProgressBar::GetInstance()->Progress();
                this->BusyCursorOff();

                //Update datastorage
                CemrgCommonUtils::AddToStorage(outputImage, imgNode->GetName(), this->GetDataStorage());
                this->GetDataStorage()->Remove(imgNode);
                mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

            }//_if
        } else
            return;
    } else
        return;
}

void EASIView::SegmentIMGS() {

    int reply = QMessageBox::question(
                NULL, "Question", "Do you have a segmentation to load?",
                QMessageBox::Yes, QMessageBox::No);

    if (reply == QMessageBox::Yes) {

        QString path = QFileDialog::getOpenFileName(
                    NULL, "Open Segmentation file", mitk::IOUtil::GetProgramPath().c_str(),
                    QmitkIOUtil::GetFileOpenFilterString());
        if (path.isEmpty())
            return;
        mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

    } else {
        //Show the plugin
        this->GetSite()->GetPage()->ShowView("org.mitk.views.segmentation");
    }//_if
}

void EASIView::BrowseM() {

    QString para = "";
    para = QFileDialog::getOpenFileName(
                NULL, "Open text file containing parameters",
                directory, QmitkIOUtil::GetFileOpenFilterString());
    m_UIMeshing.lineEdit_1->setText(para);
}

void EASIView::CreateMesh() {

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

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please select a segmentation from the Data Manager to create a mesh!");
        return;
    }

    //Find the selected node
    QString path;
    QString templatePath;
    mitk::Point3D origin;
    mitk::DataNode::Pointer segNode = nodes.at(0);
    mitk::BaseData::Pointer data = segNode->GetData();
    if (data) {

        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {

            origin = image->GetGeometry()->GetOrigin();
            int dimensions = image->GetDimension(0)*image->GetDimension(1)*image->GetDimension(2);

            try {
                //Convert image to right type
                itk::Image<uint8_t,3>::Pointer itkImage = itk::Image<uint8_t,3>::New();
                mitk::CastToItkImage(image, itkImage);
                mitk::CastToMitkImage(itkImage, image);

                //Access image volume
                mitk::ImagePixelReadAccessor<uint8_t,3> readAccess(image);
                uint8_t* pv = (uint8_t*)readAccess.GetData();

                //Prepare header of inr file (BUGS IN RELEASE MODE DUE TO NULL TERMINATOR \0)
                char header[256] = {};
                int bitlength = 8;
                const char* btype = "unsigned fixed";
                mitk::Vector3D spacing = image->GetGeometry()->GetSpacing();
                int n = sprintf(header,"#INRIMAGE-4#{\nXDIM=%d\nYDIM=%d\nZDIM=%d\nVDIM=1\nTYPE=%s\nPIXSIZE=%d bits\nCPU=decm\nVX=%6.4f\nVY=%6.4f\nVZ=%6.4f\n",image->GetDimension(0),image->GetDimension(1),image->GetDimension(2),btype,bitlength,spacing.GetElement(0),spacing.GetElement(1),spacing.GetElement(2));
                for (int i = n; i < 252; i++)
                    header[i] = '\n';

                header[252] = '#';
                header[253] = '#';
                header[254] = '}';
                header[255] = '\n';

                //Write to binary file
                std::string path = (directory + mitk::IOUtil::GetDirectorySeparator() + "converted.inr").toStdString();
                ofstream myFile(path, ios::out | ios::binary);
                myFile.write((char*)header, 256 * sizeof(char));
                myFile.write((char*)&(*pv), dimensions * sizeof(uint8_t));
                myFile.close();

                //Ask for user input to set the parameters
                QDialog* inputs = new QDialog(0,0);

                m_UIMeshing.setupUi(inputs);
                connect(m_UIMeshing.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
                connect(m_UIMeshing.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
                connect(m_UIMeshing.pushButton_1, SIGNAL(clicked()), this, SLOT(BrowseM()));

                int dialogCode = inputs->exec();

                //Act on dialog return code
                if (dialogCode == QDialog::Accepted) {

                    templatePath = m_UIMeshing.lineEdit_1->text();

                    //Checking input files
                    if (templatePath.isEmpty()) {
                        templatePath = CemrgCommonUtils::M3dlibParamFileGenerator(directory);
                    }//_if

                    //Run Mesh3DTool
                    this->BusyCursorOn();
                    mitk::ProgressBar::GetInstance()->AddStepsToDo(2);
                    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
                    cmd->SetUseDockerContainersOff(); // meshtools3d attempts the local libraries and falls back to docker.
                    cmd->SetDockerImage(QString("alonsojasl/meshtools3d:v1.0"));
                    QString output = cmd->ExecuteCreateCGALMesh(directory, "CGALMesh", templatePath);
                    QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
                    mitk::ProgressBar::GetInstance()->Progress();

                    //Prepare Visualisation
                    mitk::BaseData::Pointer meshData = mitk::IOUtil::Load(output.toStdString()).at(0);
                    mitk::UnstructuredGrid::Pointer mitkVtkGrid = dynamic_cast<mitk::UnstructuredGrid*>(meshData.GetPointer());
                    vtkSmartPointer<vtkUnstructuredGrid> vtkGrid = mitkVtkGrid->GetVtkUnstructuredGrid();
                    for (vtkIdType i=0; i<vtkGrid->GetNumberOfPoints(); i++) {
                        double* point = vtkGrid->GetPoint(i);
                        point[0] += origin.GetElement(0);
                        point[1] += origin.GetElement(1);
                        point[2] += origin.GetElement(2);
                        vtkGrid->GetPoints()->SetPoint(i, point);
                    }//_for
                    CemrgCommonUtils::AddToStorage(mitkVtkGrid, "CGALMesh", this->GetDataStorage());
                    mitk::ProgressBar::GetInstance()->Progress();
                    this->BusyCursorOff();

                } else if (dialogCode == QDialog::Rejected) {
                    inputs->close();
                    inputs->deleteLater();
                }//_if
            } catch(mitk::Exception& e) {
                //Deal with the situation not to have access
                qDebug() << e.GetDescription();
                return;
            }//_try
        } else
            return;
    } else
        return;
}

void EASIView::ActivationSites() {

    //Toggle visibility of buttons
    if (m_Controls.button_5_1->isVisible()) {
        m_Controls.button_5_1->setVisible(false);
    } else {
        m_Controls.button_5_1->setVisible(true);
        //Show the plugin
        QMessageBox::information(
                    NULL, "Attention",
                    "Please select the activation site and then press Confirm Site button.");
        this->GetSite()->GetPage()->ShowView("org.mitk.views.pointsetinteraction");
    }
}

void EASIView::ConfrmSITE() {

    //Check for selection of mesh
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please select a pointset from the Data Manager to confirm the activation site!");
        return;
    }//_if

    bool ok;
    mitk::DataNode::Pointer pstNode = nodes.at(0);
    mitk::PointSet::Pointer pstObje = dynamic_cast<mitk::PointSet*>(pstNode->GetData());
    if (pstObje && pstObje->GetSize() != 0) {

        //Ask for the size
        int scaleFactor = QInputDialog::getInt(NULL, tr("Activation Site"), tr("Enter the size:"), 5, 1, 10, 1, &ok);
        if (!ok) {
            QMessageBox::warning(NULL, "Attention", "Please enter a valid value for the size of activation site!");
            return;
        }

        //Prepare activation sphere
        mitk::DataNode::Pointer actiNode = mitk::DataNode::New();
        mitk::Ellipsoid::Pointer actiObject = mitk::Ellipsoid::New();
        actiNode->SetData(actiObject);
        actiNode->SetProperty("opacity", mitk::FloatProperty::New(0.4));
        actiNode->SetProperty("color", mitk::ColorProperty::New(1.0, 0.0, 0.0));
        actiNode->SetProperty("name", mitk::StringProperty::New("Activation Sites"));
        this->GetDataStorage()->Add(actiNode);

        //Set size of activation site
        mitk::Point3D scale;
        scale[0] = scaleFactor; scale[1] = scaleFactor; scale[2] = scaleFactor;
        actiObject->SetOrigin(pstObje->GetPoint(0));
        mitk::Point3D anchorPoint = actiObject->GetGeometry()->GetCenter();
        mitk::ScaleOperation* doOp = new mitk::ScaleOperation(mitk::OpSCALE, scale, anchorPoint);
        actiObject->GetGeometry()->ExecuteOperation(doOp);
        mitk::RenderingManager::GetInstance()->RequestUpdateAll();

    } else {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please select a pointset from the Data Manager to confirm the activation site!");
        return;
    }//_if
}

//#include <itkImageDuplicator.h>
//#include <itkLabelShapeKeepNObjectsImageFilter.h>
//#include <mitkIOUtil.h>
//#include <itkImageRegionIteratorWithIndex.h>
//#include <itkConnectedComponentImageFilter.h>
//#include <mitkImageCast.h>
//#include <mitkITKImageImport.h>
void EASIView::Simulation() {

//    typedef itk::Image<short, 3> ImageTypeCHAR;
//    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
//    for (int i=11; i<=30; i++) {

        //Prepare dirctories (1)
//        QString direct = "/home/or15/Downloads/Test";

        //Get LGE and MRA image (2)
//       QString cnnPath = direct + "/out_" + QString::number(i) + "_3D.nii.gz";

        //Clean the segmentation
//        typedef itk::ImageRegionIteratorWithIndex<ImageTypeCHAR> ItType;
//        ImageTypeCHAR::Pointer orgSegImage = ImageTypeCHAR::New();
//        CastToItkImage(mitk::IOUtil::Load<mitk::Image>(cnnPath.toStdString()), orgSegImage);
//        typedef itk::ConnectedComponentImageFilter<ImageTypeCHAR, ImageTypeCHAR> ConnectedComponentImageFilterType;
//        ConnectedComponentImageFilterType::Pointer connected1 = ConnectedComponentImageFilterType::New();
//        connected1->SetInput(orgSegImage);
//        connected1->Update();
//        typedef itk::LabelShapeKeepNObjectsImageFilter<ImageTypeCHAR> LabelShapeKeepNObjImgFilterType;
//        LabelShapeKeepNObjImgFilterType::Pointer lblShpKpNObjImgFltr1 = LabelShapeKeepNObjImgFilterType::New();
//        lblShpKpNObjImgFltr1->SetInput(connected1->GetOutput());
//        lblShpKpNObjImgFltr1->SetBackgroundValue(0);
//        lblShpKpNObjImgFltr1->SetNumberOfObjects(1);
//        lblShpKpNObjImgFltr1->SetAttribute(LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
//        lblShpKpNObjImgFltr1->Update();
//        using DuplicatorType = itk::ImageDuplicator<ImageTypeCHAR>;
//        DuplicatorType::Pointer duplicator = DuplicatorType::New();
//        duplicator->SetInputImage(lblShpKpNObjImgFltr1->GetOutput());
//        duplicator->Update();
//        ItType itDUP(duplicator->GetOutput(), duplicator->GetOutput()->GetRequestedRegion());
//        for (itDUP.GoToBegin(); !itDUP.IsAtEnd(); ++itDUP)
//            if ((int)itDUP.Get() != 0)
//                itDUP.Set(1);
//        QString segCleanPath = direct + "/prodClean_" + QString::number(i) + ".nii";
//        mitk::IOUtil::Save(mitk::ImportItkImage(duplicator->GetOutput()), segCleanPath.toStdString());
//   }
}

void EASIView::LoadMesh() {

    QString path = "";
    path = QFileDialog::getOpenFileName(
                NULL, "Open Mesh Data File",
                directory, QmitkIOUtil::GetFileOpenFilterString());
    CemrgCommonUtils::AddToStorage(
                CemrgCommonUtils::LoadVTKMesh(path.toStdString()), "Mesh", this->GetDataStorage());
}

void EASIView::Reset() {

    try {

        ctkPluginContext* context = mitk::kcl_cemrgapp_easi_Activator::getContext();
        mitk::IDataStorageService* dss = 0;
        ctkServiceReference dsRef = context->getServiceReference<mitk::IDataStorageService>();

        if (dsRef)
            dss = context->getService<mitk::IDataStorageService>(dsRef);

        if (!dss) {
            MITK_WARN << "IDataStorageService service not available. Unable to close project.";
            context->ungetService(dsRef);
            return;
        }

        mitk::IDataStorageReference::Pointer dataStorageRef = dss->GetActiveDataStorage();
        if (dataStorageRef.IsNull()) {
            //No active data storage set (i.e. not editor with a DataStorageEditorInput is active).
            dataStorageRef = dss->GetDefaultDataStorage();
        }

        mitk::DataStorage::Pointer dataStorage = dataStorageRef->GetDataStorage();
        if (dataStorage.IsNull()) {
            MITK_WARN << "No data storage available. Cannot close project.";
            return;
        }

        //Check if we got the default datastorage and if there is anything else then helper object in the storage
        if (dataStorageRef->IsDefault() && dataStorage->GetSubset(
                    mitk::NodePredicateNot::New(
                        mitk::NodePredicateProperty::New("helper object", mitk::BoolProperty::New(true))))->empty())
            return;

        //Remove everything
        mitk::DataStorage::SetOfObjects::ConstPointer nodesToRemove = dataStorage->GetAll();
        dataStorage->Remove(nodesToRemove);

        //Remove the datastorage from the data storage service
        dss->RemoveDataStorageReference(dataStorageRef);

        //Close all editors with this data storage as input
        mitk::DataStorageEditorInput::Pointer dsInput(new mitk::DataStorageEditorInput(dataStorageRef));
        QList<berry::IEditorReference::Pointer> dsEditors =
                this->GetSite()->GetPage()->FindEditors(dsInput, QString(), berry::IWorkbenchPage::MATCH_INPUT);

        if (!dsEditors.empty()) {
            QList<berry::IEditorReference::Pointer> editorsToClose = dsEditors;
            this->GetSite()->GetPage()->CloseEditors(editorsToClose, false);
        }

    } catch (std::exception& e) {

        MITK_ERROR << "Exception caught during closing project: " << e.what();
        QMessageBox::warning(NULL, "Error", QString("An error occurred during Close Project: %1").arg(e.what()));
    }//_try

    //Clear project directory
    directory.clear();
    m_Controls.button_2_2->setText("Crop Images");
}

void EASIView::LoadFibresView(){
    std::cout << "LoadFibresView" << '\n';
    if(!RequestProjectDirectoryFromUser()) return;

    //Show the plugin
    this->GetSite()->GetPage()->ResetPerspective();
    FibresView::SetDirectoryFile(directory, "CGALMesh");
    this->GetSite()->GetPage()->ShowView("org.mitk.views.fibres");
}

// helper functions
bool EASIView::RequestProjectDirectoryFromUser(){
    bool succesfulAssignment = true;
    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        MITK_INFO << "Directory is empty. Requesting user for directory.";
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        MITK_INFO << ("Directory selected:" + directory).toStdString();
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            MITK_WARN << "Please select a project directory with no spaces in the path!";
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            succesfulAssignment = false;
        }//_if
    } else {
        MITK_INFO << ("Project directory already set: " + directory).toStdString();
    }

    return succesfulAssignment;
}
