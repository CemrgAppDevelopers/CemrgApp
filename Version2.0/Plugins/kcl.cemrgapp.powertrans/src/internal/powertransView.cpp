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
 * Power Transmitter Calculations Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * angela.lee@kcl.ac.uk
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
#include <mitkStandaloneDataStorage.h>
#include <mitkImage.h>

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>
#include <QSignalMapper>

// Micro services
#include <usModuleRegistry.h>

// VTK
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

// CemrgAppModule
#include <CemrgCommandLine.h>
#include <CemrgCommonUtils.h>
#include <CemrgAhaUtils.h>
#include <CemrgPower.h>
#include "kcl_cemrgapp_powertrans_Activator.h"
#include "powertransView.h"

// Generic
#include <numeric>
#include <iostream>

const std::string powertransView::VIEW_ID = "org.mitk.views.powertrans";

void powertransView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));
    connect(m_Controls.button_3, SIGNAL(clicked()), this, SLOT(VolumeRendering()));
    connect(m_Controls.button_4, SIGNAL(clicked()), this, SLOT(MapPowerTop()));
    connect(m_Controls.button_5, SIGNAL(clicked()), this, SLOT(CalcPowerTop()));
    connect(m_Controls.button_6, SIGNAL(clicked()), this, SLOT(MapAHATop()));
    connect(m_Controls.button_r, SIGNAL(clicked()), this, SLOT(Reset()));
    //Set visibility of buttons
    m_Controls.button_2_1->setVisible(false);
    m_Controls.button_2_2->setVisible(false);
    m_Controls.button_2_3->setVisible(false);
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));
    connect(m_Controls.button_2_2, SIGNAL(clicked()), this, SLOT(CropinIMGS()));
    connect(m_Controls.button_2_3, SIGNAL(clicked()), this, SLOT(ResampIMGS()));

    m_Controls.button_4_1->setVisible(false);
    m_Controls.button_4_2->setVisible(false);
    m_Controls.button_4_3->setVisible(false);
    connect(m_Controls.button_4_1, SIGNAL(clicked()), this, SLOT(ResetRibSpacing()));
    connect(m_Controls.button_4_2, SIGNAL(clicked()), this, SLOT(LandmarkSelection()));
    connect(m_Controls.button_4_3, SIGNAL(clicked()), this, SLOT(MapPowerTransLM()));

    m_Controls.button_5_1->setVisible(false);
    m_Controls.button_5_2->setVisible(false);
    connect(m_Controls.button_5_1, SIGNAL(clicked()), this, SLOT(LoadMesh()));
    connect(m_Controls.button_5_2, SIGNAL(clicked()), this, SLOT(CalculatePower()));

    m_Controls.button_6_1->setVisible(false);
    m_Controls.button_6_2->setVisible(false);
    m_Controls.button_6_3->setVisible(false);
    connect(m_Controls.button_6_1, SIGNAL(clicked()), this, SLOT(AHALandmarkSelection()));
    connect(m_Controls.button_6_2, SIGNAL(clicked()), this, SLOT(MapAHA()));
    connect(m_Controls.button_6_3, SIGNAL(clicked()), this, SLOT(MapAHAfromInput()));

}

void powertransView::SetFocus() {

    m_Controls.button_1->setFocus();
}

void powertransView::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void powertransView::LoadDICOM() {

    //Use MITK DICOM editor
    QString editor_id = "org.mitk.editors.dicomeditor";
    berry::IEditorInput::Pointer input(new berry::FileEditorInput(QString()));
    this->GetSite()->GetPage()->OpenEditor(input, editor_id);

}

void powertransView::ProcessIMGS() {

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

void powertransView::ConvertNII() {

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
        }
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

void powertransView::CropinIMGS() {

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

void powertransView::ResampIMGS() {

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

void powertransView::VolumeRendering() {
    //Show the plugin
    this->GetSite()->GetPage()->ShowView("org.mitk.views.volumevisualization");

}

void powertransView::MapPowerTop() {

    //Toggle visibility of buttons
    if (m_Controls.button_4_1->isVisible()) {
        m_Controls.button_4_1->setVisible(false);
        m_Controls.button_4_2->setVisible(false);
        m_Controls.button_4_3->setVisible(false);
    } else {
        m_Controls.button_4_1->setVisible(true);
        m_Controls.button_4_2->setVisible(true);
        m_Controls.button_4_3->setVisible(true);
    }
}

void powertransView::ResetRibSpacing() {
    QDialog* inputs = new QDialog(0,0);
    m_UIRibSpacing.setupUi(inputs);

    connect(m_UIRibSpacing.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UIRibSpacing.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));

    int dialogCode = inputs->exec();
    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {

        QString ribSpacingStr = m_UIRibSpacing.lineEdit_1->text();

        //Checking input rib spacing
        if (ribSpacingStr.isEmpty()) {
            ribSpacing=5;
            QMessageBox::warning(NULL, "Attention", "Default: Rib spacing 5 is used");

        } else {
            bool flag;
            ribSpacing=ribSpacingStr.toInt(&flag);
            if (!flag) {
                QMessageBox::warning(NULL, "Attention", "Wrong input for rib spacing. Using default rib spacing = 5");
                ribSpacing=5;
            }
        }
        inputs->deleteLater();
    } else if (dialogCode == QDialog::Rejected) {
        inputs->close();
        inputs->deleteLater();
        return;
    }
}

void powertransView::LandmarkSelection() {
    /*
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

    //Ask for user input to set the ribSpacing number
    if (ribSpacing==0) {
        QDialog* inputs = new QDialog(0,0);
        m_UIRibSpacing.setupUi(inputs);

        connect(m_UIRibSpacing.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
        connect(m_UIRibSpacing.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));

        int dialogCode = inputs->exec();
        //Act on dialog return code
        if (dialogCode == QDialog::Accepted) {

            QString ribSpacingStr = m_UIRibSpacing.lineEdit_1->text();

            //Checking input rib spacing
            if (ribSpacingStr.isEmpty()) {
                ribSpacing=5;
                QMessageBox::warning(NULL, "Attention", "Default: Rib spacing 5 is used. \n\nPlease select 2 or 3 landmark points:\n\n1. Mid-medial edge of power transmitter\n2. 10mm along the length of the transmitter (RHS of landmark 1)\n3. [opt] 10mm along the width of the transmitter (below landmark 1)");

            } else {
                bool flag;
                ribSpacing=ribSpacingStr.toInt(&flag);
                if (!flag) {
                    QMessageBox::warning(NULL, "Attention", "Wrong input for rib spacing. Using default rib spacing = 5");
                    ribSpacing=5;
                }
            }
            inputs->deleteLater();
        } else if (dialogCode == QDialog::Rejected) {
            inputs->close();
            inputs->deleteLater();
            return;
        }
    }
    */
    QMessageBox::information(
                NULL, "Attention",
                "Please select 2 or 3 landmark points:\n\n 1. Mid-medial edge of power transmitter\n 2. 10mm along the length of the transmitter (RHS of landmark 1)\n 3. [opt] 10mm along the width of the transmitter (below landmark 1)");

    this->GetSite()->GetPage()->ShowView("org.mitk.views.pointsetinteraction");

}

void powertransView::MapPowerTransLM() {

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
    //Ask for user input to set the ribSpacing number
    if (ribSpacing==0) {
        QDialog* inputs = new QDialog(0,0);
        m_UIRibSpacing.setupUi(inputs);

        connect(m_UIRibSpacing.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
        connect(m_UIRibSpacing.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));

        int dialogCode = inputs->exec();
        //Act on dialog return code
        if (dialogCode == QDialog::Accepted) {

            QString ribSpacingStr = m_UIRibSpacing.lineEdit_1->text();

            //Checking input rib spacing
            if (ribSpacingStr.isEmpty()) {
                ribSpacing=5;
                QMessageBox::information(NULL, "Attention", "Default: Rib spacing 5 is used");

            } else {
                bool flag;
                ribSpacing=ribSpacingStr.toInt(&flag);
                if (!flag) {
                    QMessageBox::warning(NULL, "Attention", "Wrong input for rib spacing. Using default rib spacing = 5");
                    ribSpacing=5;
                }
            }
            inputs->deleteLater();
        } else if (dialogCode == QDialog::Rejected) {
            inputs->close();
            inputs->deleteLater();
            return;
        }
    }

    //Check for selection of landmarks
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if ((nodes.empty()) || (!dynamic_cast<mitk::PointSet*>(nodes.front()->GetData()))) {
        QMessageBox::warning(NULL, "Attention", "Please select landmarks from the Data Manager to continue!");
        return;
        //} else {
        //    QString path = directory + mitk::IOUtil::GetDirectorySeparator() + "Ribspace" + ribSpacing + ".csv";
        //    mitk::IOUtil::Save(lmNode, path.toStdString());
    }
    //mitk::DataNode::Pointer lmNode = nodes.front();

    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
    power = std::unique_ptr<CemrgPower>(new CemrgPower(directory, ribSpacing));
    mitk::Surface::Pointer surface = power->MapPowerTransmitterToLandmarks(nodes.front());
    if (surface->IsEmpty()) {
        mitk::ProgressBar::GetInstance()->Progress();
        this->BusyCursorOff();
        QMessageBox::warning(NULL, "Attention", "Was unable to Map Power Transmitter to given landmarks!");
        return;
    }
    CemrgCommonUtils::AddToStorage(surface, "PowerTransmitter", this->GetDataStorage());
    mitk::ProgressBar::GetInstance()->Progress();
    this->BusyCursorOff();
}

// void powertransView::ConfirmSite() {
/*
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
    //Ask for user input to set the ribSpacing number
    if (ribSpacing==0) {
        QDialog* inputs = new QDialog(0,0);
        m_UIRibSpacing.setupUi(inputs);

        connect(m_UIRibSpacing.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
        connect(m_UIRibSpacing.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));

        int dialogCode = inputs->exec();
        //Act on dialog return code
        if (dialogCode == QDialog::Accepted) {

            QString ribSpacingStr = m_UIRibSpacing.lineEdit_1->text();

            //Checking input rib spacing
            if (ribSpacingStr.isEmpty()) {
                ribSpacing=5;
                QMessageBox::information(NULL, "Attention", "Default: Rib spacing 5 is used");

            } else {
                bool flag;
                ribSpacing=ribSpacingStr.toInt(&flag);
                if (!flag) {
                    QMessageBox::warning(NULL, "Attention", "Wrong input for rib spacing. Using default rib spacing = 5");
                    ribSpacing=5;
                }
            }
            inputs->deleteLater();
        } else if (dialogCode == QDialog::Rejected) {
            inputs->close();
            inputs->deleteLater();
            return;
        }
    }

    //Check for selection of landmarks
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();

    if ((nodes.empty()) || (!dynamic_cast<mitk::PointSet*>(nodes.front()->GetData()))) {
        QMessageBox::warning(NULL, "Attention", "Please select landmarks from the Data Manager to save sites!");
        return;
    } else {
        //mitk::DataNode::Pointer lmNode = nodes.front();
        QString path = directory + mitk::IOUtil::GetDirectorySeparator() + "Ribspace" + ribSpacing + ".mps";
        //mitk::IOUtil::Save(nodes.front(), path.toStdString());
    }
*/
//}

void powertransView::CalcPowerTop() {

    //Toggle visibility of buttons
    if (m_Controls.button_5_1->isVisible()) {
        m_Controls.button_5_1->setVisible(false);
        m_Controls.button_5_2->setVisible(false);
    } else {
        m_Controls.button_5_1->setVisible(true);
        m_Controls.button_5_2->setVisible(true);
    }
}
void powertransView::LoadMesh() {

    QString path = "";
    path = QFileDialog::getOpenFileName(
                NULL, "Open Mesh Data File",
                directory, QmitkIOUtil::GetFileOpenFilterString());
    if (path.isEmpty() || path.simplified().contains(" ")) {
        QMessageBox::warning(NULL, "Attention", "Please select endo mesh file to open!");
        directory = QString();
        return;
    }//_if

    CemrgCommonUtils::AddToStorage(
                CemrgCommonUtils::LoadVTKMesh(path.toStdString()), "Mesh", this->GetDataStorage());
}

void powertransView::CalculatePower() {

    //
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
    //Ask for user input to set the ribSpacing number
    if (ribSpacing==0) {
        QDialog* inputs = new QDialog(0,0);
        m_UIRibSpacing.setupUi(inputs);

        connect(m_UIRibSpacing.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
        connect(m_UIRibSpacing.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));

        int dialogCode = inputs->exec();
        //Act on dialog return code
        if (dialogCode == QDialog::Accepted) {

            QString ribSpacingStr = m_UIRibSpacing.lineEdit_1->text();

            //Checking input rib spacing
            if (ribSpacingStr.isEmpty()) {
                ribSpacing=5;
                QMessageBox::warning(NULL, "Attention", "Default: Rib spacing 5 is used");

            } else {
                bool flag;
                ribSpacing=ribSpacingStr.toInt(&flag);
                if (!flag) {
                    QMessageBox::warning(NULL, "Attention", "Wrong input for rib spacing. Using default rib spacing = 5");
                    ribSpacing=5;
                }
            }
            inputs->deleteLater();
        } else if (dialogCode == QDialog::Rejected) {
            inputs->close();
            inputs->deleteLater();
            return;
        }
    }

    //Check for selection of endo mesh
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(NULL, "Attention", "Please select input mesh from the Data Manager to calculate power!");
        return;
    }
    mitk::DataNode::Pointer endoNode = nodes.at(0);
    mitk::BaseData::Pointer data = endoNode->GetData();
    if (data) {
        mitk::Surface::Pointer surface = dynamic_cast<mitk::Surface*>(data.GetPointer());
        if (surface) {

            QString outstr="Calculated power on selected mesh for Rib spacing" + QString::number(ribSpacing);
            QMessageBox::information(NULL, "Attention", outstr);

            this->BusyCursorOn();
            mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
            power = std::unique_ptr<CemrgPower>(new CemrgPower(directory, ribSpacing));
            mitk::Surface::Pointer outputEndoMesh = power->CalculateAcousticIntensity(surface);
            // Load in mesh
            //AddToStorage("PowerMap", outputEndoMesh);
            mitk::ProgressBar::GetInstance()->Progress();
            this->BusyCursorOff();

        } else {
            QMessageBox::warning(NULL, "Attention", "Please select input mesh from the Data Manager to calculate power!");
            return;
        }//_if

    } else {
        QMessageBox::warning(NULL, "Attention", "Can't pick up node from data manager");
        return;
    }//_if
}

void powertransView::MapAHATop() {

    //Toggle visibility of buttons
    if (m_Controls.button_6_1->isVisible()) {
        m_Controls.button_6_1->setVisible(false);
        m_Controls.button_6_2->setVisible(false);
        m_Controls.button_6_3->setVisible(false);

    } else {
        m_Controls.button_6_1->setVisible(true);
        m_Controls.button_6_2->setVisible(true);
        m_Controls.button_6_3->setVisible(true);

    }
}
void powertransView::AHALandmarkSelection() {

    QMessageBox::information(
                NULL, "Attention",
                "Please select 6 points in order:\n\n1 on the Apex\n3 on the Mitral Valve surface\n2 on the Right Ventricle cusps\n");

    this->GetSite()->GetPage()->ShowView("org.mitk.views.pointsetinteraction");

}
void powertransView::MapAHAfromInput() {

    //Ask for user input to set the parameters
    QDialog* inputs = new QDialog(0,0);
    QSignalMapper* signalMapper = new QSignalMapper(this);

    m_UIAhaInput.setupUi(inputs);
    connect(m_UIAhaInput.pushButton_1, SIGNAL(clicked()), signalMapper, SLOT(map()));
    connect(m_UIAhaInput.pushButton_2, SIGNAL(clicked()), signalMapper, SLOT(map()));
    signalMapper->setMapping(m_UIAhaInput.pushButton_1, "1"+directory + mitk::IOUtil::GetDirectorySeparator() + "ebr" + QString::number(ribSpacing) + ".vtk");
    signalMapper->setMapping(m_UIAhaInput.pushButton_2, "2"+directory );
    connect(signalMapper, SIGNAL(mapped(QString)), this, SLOT(BrowseA(const QString&)));

    int dialogCode = inputs->exec();

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {

        //Load input mesh, dofin file
        QString endomesh = m_UIAhaInput.lineEdit_1->text();
        QString AHAlm = m_UIAhaInput.lineEdit_2->text();

        mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(endomesh.toStdString());
        mitk::PointSet::Pointer pointset = mitk::IOUtil::Load<mitk::PointSet>(AHAlm.toStdString());

        //Checking input values
        if (endomesh.isEmpty() || AHAlm.isEmpty()) {
            QMessageBox::warning(NULL, "Attention", "Please select the input files to apply AHA mapping!");
            signalMapper->deleteLater();
            inputs->deleteLater();
            return;
        }

        this->BusyCursorOn();
        mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
        power = std::unique_ptr<CemrgPower>(new CemrgPower(directory, ribSpacing));
        mitk::Surface::Pointer outputAHAMesh = power->ReferenceAHA(pointset, surface);
        CemrgCommonUtils::AddToStorage(outputAHAMesh, "PowerMap", this->GetDataStorage());
        mitk::ProgressBar::GetInstance()->Progress();
        this->BusyCursorOff();

    } else if (dialogCode == QDialog::Rejected) {

        inputs->close();
        signalMapper->deleteLater();
        inputs->deleteLater();
        return;

    }//_if
}

void powertransView::MapAHA() {

    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    //Sort the two images

    if (nodes.size() != 2) {
        QMessageBox::warning(NULL, "Attention", "1) Please select  pointset and mesh for AHA mapping");
        return;
    }

    mitk::DataNode::Pointer _1st = nodes.at(0);
    mitk::DataNode::Pointer _2nd = nodes.at(1);
    mitk::BaseData::Pointer dat1 = _1st->GetData();
    mitk::BaseData::Pointer dat2 = _2nd->GetData();

    // to get the name of the file in the Data manager.
    QString name = QString::fromStdString(_2nd->GetName());
    // create path
    QString path2surface = directory + mitk::IOUtil::GetDirectorySeparator() + name;
    // load surface
    mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(path2surface.toStdString());

    if (dat1 && dat2) {
        //Test if this data items are image
        mitk::PointSet::Pointer pointset = dynamic_cast<mitk::PointSet*>(dat1.GetPointer());
        //mitk::Surface::Pointer surface = dynamic_cast<mitk::Surface*>(dat2.GetPointer());
        //if (pointset && surface) {

        if (!pointset) {
            QMessageBox::warning(NULL, "Attention", "Cannot dynamic cast pointset");
        }
        //mitk::Surface::Pointer surface = dynamic_cast<mitk::Surface*>(nodes.at(1)->GetData());

        if (!surface) {
            QMessageBox::warning(NULL, "Attention", "Cannot dynamic cast surface");
        }

        //
        if (!pointset || !surface) {
            mitk::Surface::Pointer surface = dynamic_cast<mitk::Surface*>(nodes.at(0)->GetData());
            mitk::DataNode::Pointer pointset = dynamic_cast<mitk::DataNode*>(nodes.at(1)->GetData());
            if (!pointset || !surface) {
                QMessageBox::warning(NULL, "Attention", "Cannot dynamic cast pointset and mesh for AHA mapping");
                return;
            }
        }

        this->BusyCursorOn();
        mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
        power = std::unique_ptr<CemrgPower>(new CemrgPower(directory, ribSpacing));
        mitk::Surface::Pointer outputAHAMesh = power->ReferenceAHA(pointset, surface);
        CemrgCommonUtils::AddToStorage(outputAHAMesh, "PowerMap", this->GetDataStorage());
        mitk::ProgressBar::GetInstance()->Progress();
        this->BusyCursorOff();

    } else {
        QMessageBox::warning(NULL, "Attention", "2) Please select pointset and mesh for AHA mapping");
        return;
    }//_if
}

void powertransView::Reset() {

    try {

        ctkPluginContext* context = mitk::kcl_cemrgapp_powertrans_Activator::getContext();
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
}
