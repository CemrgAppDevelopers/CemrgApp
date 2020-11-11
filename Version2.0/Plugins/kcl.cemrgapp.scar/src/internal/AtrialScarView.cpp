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
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Blueberry
#include <berryFileEditorInput.h>
#include <berryIWorkbenchPage.h>
#include <berryISelectionProvider.h>
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

// Qmitk
#include <mitkImage.h>
#include <QmitkIOUtil.h>
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkLookupTableProperty.h>
#include <mitkVtkScalarModeProperty.h>
#include <mitkIDataStorageService.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkDataStorageEditorInput.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include "kcl_cemrgapp_scar_Activator.h"
#include "AtrialScarView.h"
#include "AtrialScarClipperView.h"
#include "ScarCalculationsView.h"

// VTK
#include <vtkPolyData.h>
#include <vtkLookupTable.h>
#include <vtkSphereSource.h>
#include <vtkImageResize.h>
#include <vtkImageChangeInformation.h>
#include <vtkCellDataToPointData.h>
#include <vtkClipPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkDecimatePro.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkTimerLog.h>
#include <vtkPPolyDataNormals.h>

// ITK
#include <itkResampleImageFilter.h>
#include <itkBinaryCrossStructuringElement.h>
#include <itkGrayscaleErodeImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkLabelShapeKeepNObjectsImageFilter.h>
#include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkImageFileWriter.h>

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>
#include <QDir>
#include <QDirIterator>
#include <QFileInfo>
#include <QStringList>

// CemrgAppModule
#include <CemrgAtriaClipper.h>
#include <CemrgCommandLine.h>
#include <CemrgMeasure.h>
#include <CemrgCommonUtils.h>
#include <numeric>

const std::string AtrialScarView::VIEW_ID = "org.mitk.views.scar";

void AtrialScarView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));
    connect(m_Controls.button_3, SIGNAL(clicked()), this, SLOT(AnalysisChoice()));
    connect(m_Controls.button_4, SIGNAL(clicked()), this, SLOT(SegmentIMGS()));
    connect(m_Controls.button_5, SIGNAL(clicked()), this, SLOT(CreateSurf()));
    connect(m_Controls.button_6, SIGNAL(clicked()), this, SLOT(ScarMap()));
    connect(m_Controls.button_7, SIGNAL(clicked()), this, SLOT(Threshold()));
    connect(m_Controls.button_x, SIGNAL(clicked()), this, SLOT(Registration()));
    connect(m_Controls.button_y, SIGNAL(clicked()), this, SLOT(ClipPVeins()));
    connect(m_Controls.button_z, SIGNAL(clicked()), this, SLOT(ClipperMV()));
    connect(m_Controls.button_s, SIGNAL(clicked()), this, SLOT(Sphericity()));
    connect(m_Controls.button_c, SIGNAL(clicked()), this, SLOT(ExtraCalcs()));
    connect(m_Controls.button_r, SIGNAL(clicked()), this, SLOT(ResetMain()));

    //Sub-buttons signals
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));
    connect(m_Controls.button_x_1, SIGNAL(clicked()), this, SLOT(Register()));
    connect(m_Controls.button_x_2, SIGNAL(clicked()), this, SLOT(Transform()));
    connect(m_Controls.button_z_1, SIGNAL(clicked()), this, SLOT(SelectLandmarks()));
    connect(m_Controls.button_z_2, SIGNAL(clicked()), this, SLOT(ClipMitralValve()));
    connect(m_Controls.button_deb, SIGNAL(clicked()), this, SLOT(ScarDebug()));

    //Set visibility of buttons
    m_Controls.button_4->setVisible(false);
    m_Controls.button_5->setVisible(false);
    m_Controls.button_6->setVisible(false);
    m_Controls.button_7->setVisible(false);
    m_Controls.button_x->setVisible(false);
    m_Controls.button_y->setVisible(false);
    m_Controls.button_z->setVisible(false);
    m_Controls.button_s->setVisible(false);

    //Set visibility of sub-buttons
    m_Controls.button_2_1->setVisible(false);
    m_Controls.button_x_1->setVisible(false);
    m_Controls.button_x_2->setVisible(false);
    m_Controls.button_z_1->setVisible(false);
    m_Controls.button_z_2->setVisible(false);
    m_Controls.button_deb->setVisible(false);
}

void AtrialScarView::SetFocus() {

    m_Controls.button_1->setFocus();
}

void AtrialScarView::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void AtrialScarView::LoadDICOM() {

    int reply1 = QMessageBox::No;
#if defined(__APPLE__)
    MITK_INFO << "Ask user about alternative DICOM reader";
    reply1 = QMessageBox::question(NULL, "Question",
                                   "Use alternative DICOM reader?", QMessageBox::Yes, QMessageBox::No);
#endif

    if (reply1 == QMessageBox::Yes) {

        QString dicomFolder = QFileDialog::getExistingDirectory(NULL, "Open folder with DICOMs.", mitk::IOUtil::GetProgramPath().c_str(), QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        QString tmpNiftiFolder = cmd->DockerDicom2Nifti(dicomFolder);

        if (tmpNiftiFolder.compare("ERROR_IN_PROCESSING") != 0) {

            // add results in NIIs folder to Data Manager
            MITK_INFO << ("Conversion succesful. Intermediate NII folder: " + tmpNiftiFolder).toStdString();
            QMessageBox::information(NULL, "Information", "Conversion successful, press the Process Images button to continue.");
            QDir niftiFolder(tmpNiftiFolder);
            QStringList niftiFiles = niftiFolder.entryList();

            if (niftiFiles.size()>0) {

                QString thisFile, path;
                for(int ix=0; ix<niftiFiles.size(); ix++) {

                    // load here files
                    thisFile = niftiFiles.at(ix);
                    if (thisFile.contains(".nii", Qt::CaseSensitive)) {
                        if (thisFile.contains("lge", Qt::CaseInsensitive) ||  thisFile.contains("mra", Qt::CaseInsensitive)) {

                            path = niftiFolder.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + thisFile;
                            mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(path.toStdString());
                            std::string key = "dicom.series.SeriesDescription";
                            mitk::DataStorage::SetOfObjects::Pointer set = mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
                            set->Begin().Value()->GetData()->GetPropertyList()->SetStringProperty(key.c_str(), thisFile.left(thisFile.length()-4).toStdString().c_str());

                        }//_if
                    }//_if
                }//_for

            } else {

                MITK_WARN << "Problem with conversion.";
                QMessageBox::warning(NULL, "Attention", "Problem with alternative conversion. Try MITK Dicom editor?");
                return;

            }//_if
        }//_if

    } else {

        MITK_INFO << "Using MITK DICOM editor";
        QString editor_id = "org.mitk.editors.dicomeditor";
        berry::IEditorInput::Pointer input(new berry::FileEditorInput(QString()));
        this->GetSite()->GetPage()->OpenEditor(input, editor_id);

    }//_if
}

void AtrialScarView::ProcessIMGS() {

    //Toggle visibility of buttons
    if (m_Controls.button_2_1->isVisible()){
        m_Controls.button_2_1->setVisible(false);
    }
    else{
        m_Controls.button_2_1->setVisible(true);
    }
}

void AtrialScarView::ConvertNII() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 2) {
        QMessageBox::warning(NULL, "Attention", "Please load and select both LGE and CEMRA images from the Data Manager to convert!");
        return;
    }//_if

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    //Order dicoms based on their type
    std::vector<int> indexNodes;
    std::vector<std::string> seriesDscrps;
    foreach (mitk::DataNode::Pointer node, nodes) {

        std::string seriesDescription;
        node->GetData()->GetPropertyList()->GetStringProperty("dicom.series.SeriesDescription", seriesDescription);

        if (seriesDescription.find("LGE")      != seriesDescription.npos) indexNodes.push_back(0);
        else if (seriesDescription.find("MRA") != seriesDescription.npos) indexNodes.push_back(1);

        //Trim whitespaces
        seriesDescription = QString::fromStdString(seriesDescription).replace(")","").toStdString();
        seriesDescription = QString::fromStdString(seriesDescription).replace("(","").toStdString();
        seriesDescription = QString::fromStdString(seriesDescription).simplified().replace(" ","").toStdString();
        seriesDscrps.push_back(seriesDescription);
    }//_for

    //Sort indexes based on comparing values
    std::vector<int> index(indexNodes.size());
    std::iota(index.begin(), index.end(), 0);
    std::sort(index.begin(), index.end(), [&](int i1, int i2) {return indexNodes[i1]<indexNodes[i2];});

    //Warning for cases when type is not found
    size_t length1 = nodes.size();
    size_t length2 = indexNodes.size();
    bool test = std::adjacent_find(indexNodes.begin(), indexNodes.end(), std::not_equal_to<int>()) == indexNodes.end();
    if (length1 != length2 || test) {
        QMessageBox::warning(NULL, "Attention",
                             "Cannot find the type of images automatically. Revert to user order and selections in the data manager: LGE at the top, then CEMRA at the bottom!");
        index.resize(nodes.size());
        std::iota(index.begin(), index.end(), 0);
    }//_if

    //Convert to Nifti
    int ctr = 0;
    QString path, type;
    bool successfulNitfi, resampleImage, reorientToRAI;
    resampleImage = true;
    reorientToRAI = true;

    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(index.size());
    foreach (int idx, index) {
        type = (ctr==0) ? "LGE":"MRA";
        path = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + type + "-" + seriesDscrps.at(idx).c_str() + ".nii";
        successfulNitfi = CemrgCommonUtils::ConvertToNifti(nodes.at(idx)->GetData(), path, resampleImage, reorientToRAI);
        if (successfulNitfi) {
            this->GetDataStorage()->Remove(nodes.at(idx));
            std::string key = "dicom.series.SeriesDescription";
            mitk::DataStorage::SetOfObjects::Pointer set = mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
            set->Begin().Value()->GetData()->GetPropertyList()->SetStringProperty(key.c_str(), seriesDscrps.at(idx).c_str());
            ctr++;
        } else {
            mitk::ProgressBar::GetInstance()->Progress(index.size());
            return;
        }//_if
        mitk::ProgressBar::GetInstance()->Progress();
    }//for
    nodes.clear();
    this->BusyCursorOff();

    MITK_INFO << "Loading all items";
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
}

void AtrialScarView::AnalysisChoice() {

    int reply = QMessageBox::question(
                NULL, "Question", "Do you want an automatic analysis?", QMessageBox::Yes, QMessageBox::No);

    if (reply == QMessageBox::Yes) {

        MITK_INFO << "Setting up automatic analysis.";
        AutomaticAnalysis();

    } else {

        MITK_INFO << "Setting up manual analysis.";

        //Set visibility of buttons
        m_Controls.button_4->setVisible(true);
        m_Controls.button_5->setVisible(true);
        m_Controls.button_6->setVisible(true);
        m_Controls.button_7->setVisible(true);
        m_Controls.button_x->setVisible(true);
        m_Controls.button_y->setVisible(true);
        m_Controls.button_z->setVisible(true);
        m_Controls.button_s->setVisible(true);
    }//_if
}

void AtrialScarView::AutomaticAnalysis() {

    MITK_INFO << "Performing automatic analysis.";
    MITK_INFO << "============= Automatic segmentation module ====================";

    QString direct, mraPath, lgePath, cnnPath;
    bool debugging = true;

    if (directory.isEmpty()) {
        direct = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory",
                    mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        directory = direct;
    } else {
        direct = directory;
    }//_dir

    QDirIterator searchit(direct, QDirIterator::Subdirectories);

    MITK_INFO(debugging) << "[DEBUG] Searching for CEMRGNET output";

    while(searchit.hasNext()) {
        QFileInfo searchfinfo(searchit.next());
        if (searchfinfo.fileName().contains(".nii", Qt::CaseSensitive)) {
            if (searchfinfo.fileName().contains("dcm-LGE", Qt::CaseSensitive))
                lgePath = searchfinfo.absoluteFilePath();

            if (searchfinfo.fileName().contains("dcm-MRA", Qt::CaseSensitive))
                mraPath = searchfinfo.absoluteFilePath();

            if (debugging && searchfinfo.fileName().contains("LA-cemrgnet.nii", Qt::CaseSensitive))
                cnnPath = searchfinfo.absoluteFilePath();
        }
    }//_while

    QDialog* inputs = new QDialog(0,0);
    m_UIcemrgnet.setupUi(inputs);
    connect(m_UIcemrgnet.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UIcemrgnet.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    int dialogCode = inputs->exec();

    QString meType_UI;
    int minStep_UI =-1, maxStep_UI = 3;
    int methodType_UI = 2, thresh_methodType_UI = 1;
    QStringList separated_thresh_list;
    std::vector<double> values_vector;

    if (dialogCode == QDialog::Accepted) {

        MITK_INFO << "[UI] User inputs being selected.";
        MITK_INFO << "[UI] Intensity projection";
        bool ok1, ok2;
        minStep_UI = m_UIcemrgnet.minStep_lineEdit->text().toInt(&ok1);
        maxStep_UI = m_UIcemrgnet.maxStep_lineEdit->text().toInt(&ok2);
        if (!ok1) minStep_UI = -1;
        if (!ok2) maxStep_UI = 3;

        methodType_UI = m_UIcemrgnet.maxProjection_radioButton->isChecked() ? 2 : 1;
        meType_UI = m_UIcemrgnet.maxProjection_radioButton->isChecked() ? "Max" : "Mean";

        MITK_INFO << ("[UI] Using: " + meType_UI + " Intensity projection.").toStdString();
        MITK_INFO << ("[UI] In/out values: (" + QString::number(minStep_UI) + ", " +
                      QString::number(maxStep_UI) + ")").toStdString();
        MITK_INFO << QString::number(methodType_UI);
        MITK_INFO << "[UI] Thresholding information.";

        //bool ok3;
        thresh_methodType_UI = m_UIcemrgnet.iir_radioButton->isChecked() ? 1 : 2;
        QString thresh_list, whichThresh;

        if (m_UIcemrgnet.iir_radioButton->isChecked()) { // IIR method
            whichThresh = "IIR";
            thresh_list = m_UIcemrgnet.iir_textEdit->toPlainText();
            separated_thresh_list << "0.97" << "1.16";
        } else if (m_UIcemrgnet.meanSD_radioButton->isChecked()) { // SDev method
            whichThresh = "MEAN+SD";
            thresh_list = m_UIcemrgnet.meanSD_textEdit->toPlainText();
            separated_thresh_list << "2.3" << "3.3";
        }

        MITK_INFO << ("[UI] Threshold: " + whichThresh).toStdString();
        MITK_INFO << ("[UI] Threshold list: " + thresh_list).toStdString();
        MITK_INFO << QString::number(thresh_methodType_UI);

        thresh_list.remove(" ", Qt::CaseSensitive);
        if (!thresh_list.isEmpty()) {

            MITK_INFO << "[UI] Creating list of thresholds";
            separated_thresh_list.removeLast();
            separated_thresh_list.removeLast();
            separated_thresh_list = thresh_list.split("," , QString::SkipEmptyParts);
            int listspaces = separated_thresh_list.removeAll(" ");
            int listduplicates = separated_thresh_list.removeDuplicates();
            separated_thresh_list.sort();
            if (debugging) {
                MITK_INFO << ("[UI][DEBUG] Spaces: " + QString::number(listspaces)).toStdString();
                MITK_INFO << ("[UI][DEBUG] Duplicates: " + QString::number(listduplicates)).toStdString();
            }

        }//_if

        double tryNumber;
        bool vOK;
        for(int ix=0; ix<separated_thresh_list.size(); ix++) {
            MITK_INFO << separated_thresh_list.at(ix);
            tryNumber = separated_thresh_list.at(ix).toDouble(&vOK);
            if (vOK) values_vector.push_back(tryNumber);
        }

        inputs->close();
        inputs->deleteLater();

    } else if (dialogCode == QDialog::Rejected) {

        MITK_INFO << "[ATTENTION] Cancelled automatic analysis.";
        QMessageBox::warning(
                    NULL, "Automatic analysis cancelled",
                    "'Cancel' button pressed, no calculations were made.");
        inputs->close();
        inputs->deleteLater();
        return;

    }//_if

    MITK_INFO << ("Files to be read: \n\n [LGE]: " + lgePath + "\n [MRA]: " + mraPath).toStdString();

    if (!mraPath.isEmpty()) {

        vtkSmartPointer<vtkTimerLog> timerLog = vtkSmartPointer<vtkTimerLog>::New();
        typedef itk::Image<short, 3> ImageTypeSHRT;
        typedef itk::Image<short, 3> ImageTypeCHAR;
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        MITK_INFO << "[AUTOMATIC_ANALYSIS] Setting Docker on MIRTK to OFF";
        cmd->SetUseDockerContainers(_useDockerInPlugin);

        timerLog->StartTimer();
        if (cnnPath.isEmpty()) {
            cnnPath = cmd->DockerCemrgNetPrediction(mraPath);
        }

        if (!cnnPath.isEmpty()) {

            MITK_INFO << ("Successful prediction with file "+cnnPath).toStdString();
            // QString direct = finfo.absolutePath();
            MITK_INFO << "[AUTOMATIC_ANALYSIS][1] Adjust CNN label to MRA";
            mitk::Image::Pointer mraIMG = mitk::IOUtil::Load<mitk::Image>(mraPath.toStdString());
            mitk::Image::Pointer cnnIMG = mitk::IOUtil::Load<mitk::Image>(cnnPath.toStdString());
            double origin[3]; double spacing[3];
            mraIMG->GetGeometry()->GetOrigin().ToArray(origin);
            mraIMG->GetGeometry()->GetSpacing().ToArray(spacing);

            vtkSmartPointer<vtkImageResize> resizeFilter = vtkSmartPointer<vtkImageResize>::New();
            resizeFilter->SetResizeMethodToOutputDimensions();
            resizeFilter->SetOutputDimensions(mraIMG->GetDimension(0), mraIMG->GetDimension(1), mraIMG->GetDimension(2));
            resizeFilter->InterpolateOff();
            resizeFilter->SetInputData(cnnIMG->GetVtkImageData());
            resizeFilter->Update();

            vtkSmartPointer<vtkImageChangeInformation> changeFilter = vtkSmartPointer<vtkImageChangeInformation>::New();
            changeFilter->SetInputConnection(resizeFilter->GetOutputPort());
            changeFilter->SetOutputSpacing(spacing);
            changeFilter->SetOutputOrigin(origin);
            changeFilter->Update();

            mitk::Image::Pointer cnnLA = mitk::Image::New();
            cnnIMG->Initialize(changeFilter->GetOutput());
            cnnIMG->SetVolume(changeFilter->GetOutput()->GetScalarPointer());

            MITK_INFO << "[AUTOMATIC_ANALYSIS][2] Image registration";
            QString cnnPath = direct + mitk::IOUtil::GetDirectorySeparator() + "LA.nii";
            QString laregPath = direct + mitk::IOUtil::GetDirectorySeparator() + "LA-reg.nii";

            mitk::IOUtil::Save(cnnIMG, cnnPath.toStdString());
            cmd->ExecuteRegistration(direct, lgePath, mraPath); // rigid.dof is the default name
            cmd->ExecuteTransformation(direct, cnnPath, laregPath);

            MITK_INFO << "[AUTOMATIC_ANALYSIS][3] Clean segmentation";
            typedef itk::ImageRegionIteratorWithIndex<ImageTypeCHAR> ItType;
            typedef itk::ConnectedComponentImageFilter<ImageTypeCHAR, ImageTypeCHAR> ConnectedComponentImageFilterType;
            typedef itk::LabelShapeKeepNObjectsImageFilter<ImageTypeCHAR> LabelShapeKeepNObjImgFilterType;
            using DuplicatorType = itk::ImageDuplicator<ImageTypeCHAR>;

            ImageTypeCHAR::Pointer orgSegImage = ImageTypeCHAR::New();
            mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(laregPath.toStdString()), orgSegImage);

            ConnectedComponentImageFilterType::Pointer connected1 = ConnectedComponentImageFilterType::New();
            connected1->SetInput(orgSegImage);
            connected1->Update();

            LabelShapeKeepNObjImgFilterType::Pointer lblShpKpNObjImgFltr1 = LabelShapeKeepNObjImgFilterType::New();
            lblShpKpNObjImgFltr1->SetInput(connected1->GetOutput());
            lblShpKpNObjImgFltr1->SetBackgroundValue(0);
            lblShpKpNObjImgFltr1->SetNumberOfObjects(1);
            lblShpKpNObjImgFltr1->SetAttribute(
                        LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
            lblShpKpNObjImgFltr1->Update();

            DuplicatorType::Pointer duplicator = DuplicatorType::New();
            duplicator->SetInputImage(lblShpKpNObjImgFltr1->GetOutput());
            duplicator->Update();
            ItType itDUP(duplicator->GetOutput(), duplicator->GetOutput()->GetRequestedRegion());
            for (itDUP.GoToBegin(); !itDUP.IsAtEnd(); ++itDUP)
                if ((int)itDUP.Get() != 0)
                    itDUP.Set(1);
            QString segCleanPath = direct + mitk::IOUtil::GetDirectorySeparator() +
                    "prodClean.nii";
            mitk::IOUtil::Save(mitk::ImportItkImage(duplicator->GetOutput()), segCleanPath.toStdString());
            MITK_INFO << ("[...][3.1] Saved file: "+segCleanPath).toStdString();

            MITK_INFO << "[AUTOMATIC_ANALYSIS][4] Vein clipping mesh";
            QString output1 = cmd->ExecuteSurf(direct, segCleanPath, "close", 1, .5, 0, 10);
            mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(output1.toStdString());
            vtkSmartPointer<vtkDecimatePro> deci = vtkSmartPointer<vtkDecimatePro>::New();
            deci->SetInputData(shell->GetVtkPolyData());
            deci->SetTargetReduction(0.1);
            deci->PreserveTopologyOn();
            deci->Update();
            shell->SetVtkPolyData(deci->GetOutput());

            vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
            vtkSmartPointer<vtkPolyData> pd = shell->Clone()->GetVtkPolyData();
            for (int i=0; i<pd->GetNumberOfPoints(); i++) {
                double* point = pd->GetPoint(i);
                point[0] = -point[0];
                point[1] = -point[1];
                pd->GetPoints()->SetPoint(i, point);
            }//_for
            pointLocator->SetDataSet(pd);
            pointLocator->BuildLocator();

            MITK_INFO << "[AUTOMATIC_ANALYSIS][5] Separate veins";
            typedef itk::BinaryCrossStructuringElement<ImageTypeCHAR::PixelType, 3> CrossType;
            typedef itk::BinaryMorphologicalOpeningImageFilter<ImageTypeCHAR, ImageTypeCHAR, CrossType> MorphFilterType;
            typedef itk::RelabelComponentImageFilter<ImageTypeCHAR, ImageTypeCHAR> RelabelFilterType;

            ImageTypeCHAR::Pointer veinsSegImage = ImageTypeCHAR::New();
            veinsSegImage = lblShpKpNObjImgFltr1->GetOutput();
            ItType itORG(orgSegImage, orgSegImage->GetRequestedRegion());
            ItType itVEN(veinsSegImage, veinsSegImage->GetRequestedRegion());
            itORG.GoToBegin();

            for (itVEN.GoToBegin(); !itVEN.IsAtEnd(); ++itVEN) {
                if ((int)itVEN.Get() != 0)
                    itVEN.Set((int)itORG.Get());
                ++itORG;
            }
            for (itVEN.GoToBegin(); !itVEN.IsAtEnd(); ++itVEN)
                if ((int)itVEN.Get() != 2)
                    itVEN.Set(0);

            CrossType binaryCross;
            binaryCross.SetRadius(2.0);
            binaryCross.CreateStructuringElement();
            MorphFilterType::Pointer morphFilter = MorphFilterType::New();
            morphFilter->SetInput(veinsSegImage);
            morphFilter->SetKernel(binaryCross);
            morphFilter->SetForegroundValue(2);
            morphFilter->SetBackgroundValue(0);
            morphFilter->UpdateLargestPossibleRegion();
            veinsSegImage = morphFilter->GetOutput();
            mitk::IOUtil::Save(mitk::ImportItkImage(veinsSegImage), (direct + "/prodVeins.nii").toStdString());

            ConnectedComponentImageFilterType::Pointer connected2 = ConnectedComponentImageFilterType::New();
            connected2->SetInput(veinsSegImage);
            connected2->Update();

            RelabelFilterType::Pointer relabeler = RelabelFilterType::New();
            relabeler->SetInput(connected2->GetOutput());
            relabeler->Update();
            mitk::IOUtil::Save(mitk::ImportItkImage(relabeler->GetOutput()), (direct + "/prodSeparatedVeins.nii").toStdString());
            MITK_INFO << ("[...][5.1] Saved file: "+direct + "/prodSeparatedVeins.nii").toStdString();

            MITK_INFO << "[AUTOMATIC_ANALYSIS][6] Find vein landmark";
            veinsSegImage = relabeler->GetOutput();
            ItType itLMK(veinsSegImage, veinsSegImage->GetRequestedRegion());
            vtkSmartPointer<vtkIdList> pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
            pickedSeedIds->Initialize();
            std::vector<std::vector<double>> veinsCentre;
            const int nveins = static_cast<int>(connected2->GetObjectCount());

            MITK_INFO << ("[...][6.1] Number of veins found: "+QString::number(nveins)).toStdString();
            for (int j=0; j<nveins; j++) {
                int ctrVeinsVoxels = 0;
                std::vector<double> veinLandmark(3, 0.0);
                for (itLMK.GoToBegin(); !itLMK.IsAtEnd(); ++itLMK) {
                    if ((int)itLMK.Get() == (j+1)) {
                        ImageTypeCHAR::PointType point;
                        veinsSegImage->TransformIndexToPhysicalPoint(itLMK.GetIndex(), point);
                        veinLandmark[0] += point[0];
                        veinLandmark[1] += point[1];
                        veinLandmark[2] += point[2];
                        ctrVeinsVoxels++;
                    }
                }//_for
                veinLandmark[0] /= ctrVeinsVoxels;
                veinLandmark[1] /= ctrVeinsVoxels;
                veinLandmark[2] /= ctrVeinsVoxels;
                veinsCentre.push_back(veinLandmark);
            }//_nveins
            for (int j=0; j<nveins; j++) {
                double veinLandmark[3];
                veinLandmark[0] = veinsCentre.at(j)[0];
                veinLandmark[1] = veinsCentre.at(j)[1];
                veinLandmark[2] = veinsCentre.at(j)[2];
                vtkIdType id = pointLocator->FindClosestPoint(veinLandmark);
                pickedSeedIds->InsertNextId(id);
            }//_nveins
            std::vector<int> pickedSeedLabels;
            for (int j=0; j<nveins; j++)
                pickedSeedLabels.push_back(21);

            MITK_INFO << "[AUTOMATIC_ANALYSIS][7] Clip the veins";

            std::unique_ptr<CemrgAtriaClipper> clipper(new CemrgAtriaClipper(direct, shell));
            bool successful = clipper->ComputeCtrLines(pickedSeedLabels, pickedSeedIds, false);
            if (!successful) {
                QMessageBox::critical(NULL, "Attention", "Computation of Centrelines Failed!");
                return;
            }//_Check for failure
            MITK_INFO << "[...][7.1] ComputeCtrLines finished .";

            successful = clipper->ComputeCtrLinesClippers(pickedSeedLabels);
            if (!successful) {
                QMessageBox::critical(NULL, "Attention", "Computation of Clipper Planes Failed!");
                return;
            }//_if
            MITK_INFO << "[...][7.2] ComputeCtrLinesClippers finished .";

            clipper->ClipVeinsImage(pickedSeedLabels, mitk::ImportItkImage(duplicator->GetOutput()), false);
            MITK_INFO << "[...][7.3] ClipVeinsImage finished .";

            MITK_INFO << "[AUTOMATIC_ANALYSIS][8] Create a mesh from clipped segmentation of veins";
            QString output2 = cmd->ExecuteSurf(direct, (direct + "/PVeinsCroppedImage.nii"), "close", 1, .5, 0, 10);
            mitk::Surface::Pointer LAShell = mitk::IOUtil::Load<mitk::Surface>(output2.toStdString());

            MITK_INFO << "[AUTOMATIC_ANALYSIS][9] Clip the mitral valve";
            ImageTypeCHAR::Pointer mvImage = ImageTypeCHAR::New();
            mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(segCleanPath.toStdString()), mvImage);
            ItType itMVI1(mvImage, mvImage->GetRequestedRegion());
            itORG.GoToBegin();
            for (itMVI1.GoToBegin(); !itMVI1.IsAtEnd(); ++itMVI1) {
                if ((int)itMVI1.Get() != 0)
                    itMVI1.Set((int)itORG.Get());
                ++itORG;
            }
            for (itMVI1.GoToBegin(); !itMVI1.IsAtEnd(); ++itMVI1)
                if ((int)itMVI1.Get() != 3)
                    itMVI1.Set(0);
            typedef itk::ConnectedComponentImageFilter<ImageTypeCHAR, ImageTypeCHAR> ConnectedComponentImageFilterType;
            ConnectedComponentImageFilterType::Pointer connected3 = ConnectedComponentImageFilterType::New();
            connected3->SetInput(mvImage);
            connected3->Update();
            typedef itk::LabelShapeKeepNObjectsImageFilter<ImageTypeCHAR> LabelShapeKeepNObjImgFilterType;
            LabelShapeKeepNObjImgFilterType::Pointer lblShpKpNObjImgFltr2 = LabelShapeKeepNObjImgFilterType::New();
            lblShpKpNObjImgFltr2->SetInput(connected3->GetOutput());
            lblShpKpNObjImgFltr2->SetBackgroundValue(0);
            lblShpKpNObjImgFltr2->SetNumberOfObjects(1);
            lblShpKpNObjImgFltr2->SetAttribute(LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
            lblShpKpNObjImgFltr2->Update();
            mvImage = lblShpKpNObjImgFltr2->GetOutput();
            mitk::IOUtil::Save(mitk::ImportItkImage(mvImage), (direct + "/prodMVI.nii").toStdString());

            // Make vtk of prodMVI
            QString mviShellPath = cmd->ExecuteSurf(direct, "prodMVI.nii", "close", 1, 0.5, 0, 10);
            // Implement code from command line tool
            mitk::Surface::Pointer ClipperSurface = mitk::IOUtil::Load<mitk::Surface>(mviShellPath.toStdString());
            vtkSmartPointer<vtkImplicitPolyDataDistance> implicitFn = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
            implicitFn->SetInput(ClipperSurface->GetVtkPolyData());
            vtkMTimeType mtime = implicitFn->GetMTime();
            std::cout << "MTime: " << mtime<< std::endl ;
            vtkSmartPointer<vtkClipPolyData> mvclipper = vtkSmartPointer<vtkClipPolyData>::New();
            mvclipper->SetClipFunction(implicitFn);
            mvclipper->SetInputData(LAShell->GetVtkPolyData());
            mvclipper->InsideOutOff();
            mvclipper->Update();

            MITK_INFO << "[...][9.1] Extract and clean surface mesh.";
            vtkSmartPointer<vtkDataSetSurfaceFilter> surfer = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
            surfer->SetInputData(mvclipper->GetOutput());
            surfer->Update();

            MITK_INFO << "[...][9.2] Cleaning...";
            vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
            clean->SetInputConnection(surfer->GetOutputPort());
            clean->Update();

            MITK_INFO << "[...][9.3] Largest region...";
            vtkSmartPointer<vtkPolyDataConnectivityFilter> lrgRegion = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
            lrgRegion->SetInputConnection(clean->GetOutputPort());
            lrgRegion->SetExtractionModeToLargestRegion();
            lrgRegion->Update();
            clean = vtkSmartPointer<vtkCleanPolyData>::New();
            clean->SetInputConnection(lrgRegion->GetOutputPort());
            clean->Update();

            MITK_INFO << ("[...][9.4] Saving to file: " + output2).toStdString();
            LAShell->SetVtkPolyData(clean->GetOutput());
            mitk::IOUtil::Save(LAShell, output2.toStdString());

            MITK_INFO << "[AUTOMATIC_ANALYSIS][10] Scar projection";
            int minStep = minStep_UI;
            int maxStep = maxStep_UI;
            int methodType = methodType_UI;
            std::unique_ptr<CemrgScar3D> scar(new CemrgScar3D());
            scar->SetMinStep(minStep);
            scar->SetMaxStep(maxStep);
            scar->SetMethodType(methodType);
            ImageTypeCHAR::Pointer segITK = ImageTypeCHAR::New();
            mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>((direct + "/PVeinsCroppedImage.nii").toStdString()), segITK);
            ImageTypeSHRT::Pointer lgeITK = ImageTypeSHRT::New();
            mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(lgePath.toStdString()), lgeITK);
            itk::ResampleImageFilter<ImageTypeCHAR, ImageTypeCHAR>::Pointer resampleFilter;
            resampleFilter = itk::ResampleImageFilter<ImageTypeCHAR, ImageTypeCHAR>::New();
            resampleFilter->SetInput(segITK);
            resampleFilter->SetReferenceImage(lgeITK);
            resampleFilter->SetUseReferenceImage(true);
            resampleFilter->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<ImageTypeCHAR>::New());
            resampleFilter->SetDefaultPixelValue(0);
            resampleFilter->UpdateLargestPossibleRegion();
            segITK = resampleFilter->GetOutput();
            mitk::IOUtil::Save(mitk::ImportItkImage(segITK), (direct + "/PVeinsCroppedImage.nii").toStdString());
            scar->SetScarSegImage(mitk::ImportItkImage(segITK));
            mitk::Surface::Pointer scarShell = scar->Scar3D(direct.toStdString(), mitk::ImportItkImage(lgeITK));
            MITK_INFO << "[...][10.1] Converting cell to point data";
            vtkSmartPointer<vtkCellDataToPointData> cell_to_point = vtkSmartPointer<vtkCellDataToPointData>::New();
            cell_to_point->SetInputData(scarShell->GetVtkPolyData());
            cell_to_point->PassCellDataOn();
            cell_to_point->Update();
            scarShell->SetVtkPolyData(cell_to_point->GetPolyDataOutput());
            mitk::IOUtil::Save(scarShell, (direct + "/MaxScar.vtk").toStdString());
            scar->SaveScarDebugImage("Max_debugScar.nii", direct);

            MITK_INFO << "[AUTOMATIC_ANALYSIS][11] Thresholding";
            int vxls = 3;
            int threshType = thresh_methodType_UI;

            typedef itk::Image<float, 3> ImageType;
            typedef itk::BinaryBallStructuringElement<ImageTypeCHAR::PixelType, 3> BallType;
            typedef itk::GrayscaleErodeImageFilter<ImageTypeCHAR, ImageType, BallType> ErosionFilterType;
            BallType binaryBall;
            binaryBall.SetRadius(vxls);
            binaryBall.CreateStructuringElement();
            ErosionFilterType::Pointer erosionFilter = ErosionFilterType::New();
            erosionFilter->SetInput(segITK);
            erosionFilter->SetKernel(binaryBall);
            erosionFilter->UpdateLargestPossibleRegion();
            mitk::Image::Pointer roiImage = mitk::Image::New();
            roiImage = mitk::ImportItkImage(erosionFilter->GetOutput())->Clone();
            ImageType::Pointer lgeFloat = ImageType::New();
            mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(lgePath.toStdString()), lgeFloat);
            double mean = 0.0, stdv = 0.0;
            scar->CalculateMeanStd(mitk::ImportItkImage(lgeFloat), roiImage, mean, stdv);
            MITK_INFO << "[...][11.1] Creating Scar map normalised by Mean blood pool.";
            QString prodPath = direct + mitk::IOUtil::GetDirectorySeparator();
            scar->saveNormalisedScalars(mean, scarShell, (prodPath + "MaxScar_Normalised.vtk"));
            MITK_INFO << "[...][11.2] Saving to files.";
            double thisThresh, thisPercentage, thisValue;
            ofstream prodFile1, prodFileExplanation;
            prodFile1.open((prodPath + "prodThresholds.txt").toStdString());
            for(int ix=0; (unsigned) ix < values_vector.size(); ix++) {
                thisValue = values_vector.at(ix);
                thisThresh = (threshType == 1) ? mean*thisValue : mean + thisValue*stdv;
                thisPercentage = scar->Thresholding(thisThresh);
                prodFile1 << thisValue << "\n";
                prodFile1 << threshType << "\n";
                prodFile1 << mean << "\n";
                prodFile1 << stdv << "\n";
                prodFile1 << thisThresh << "\n";
                prodFile1 << "SCORE: " << thisPercentage << "\n";
                prodFile1 << "=============== separation ================\n";
            }
            prodFileExplanation.open((prodPath + "prodThresholds_Guide.txt").toStdString());
            prodFileExplanation << "VALUE\n";
            prodFileExplanation << "THRESHOLD TYPE: (1 = V*IIR, 2 = MEAN + V*STDev)\n";
            prodFileExplanation << "MEAN INTENSITY\n";
            prodFileExplanation << "STANDARD DEVIATION (STDev)\n";
            prodFileExplanation << "THRESHOLD\n";
            prodFileExplanation << "SCAR SCORE (percentage)\n";
            prodFileExplanation << "=============== separation ================";
            prodFile1.close();
            prodFileExplanation.close();
            timerLog->StopTimer();

            QStringList rtminsec = QString::number(timerLog->GetElapsedTime()/60).split(".");
            QString rtmin = rtminsec.at(0);
            QString rtsec = QString::number(("0."+rtminsec.at(1)).toFloat()*60, 'f',1);
            QString outstr = "Operation finshed in " + rtmin + " min and " + rtsec + " s.";
            MITK_INFO << "[AUTOMATIC_ANALYSIS][FINISHED]";
            QMessageBox::information(NULL, "Automatic analysis", outstr);

        } else
            QMessageBox::warning(NULL, "Attention", "Error with automatic segmentation! Check the LOG file.");
    } else
        QMessageBox::information(NULL, "Attention", "Operation Cancelled!");
}

void AtrialScarView::SegmentIMGS() {

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    int reply1 = QMessageBox::question(
                NULL, "Question", "Do you have a segmentation to load?", QMessageBox::Yes, QMessageBox::No);

    if (reply1 == QMessageBox::Yes) {
        QString path = QFileDialog::getOpenFileName(NULL, "Open Segmentation file",
                                                    directory.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());
        if (path.isEmpty()) return;
        mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

        //Restore image name
        QFileInfo fullPathInfo(path);
        fileName = fullPathInfo.fileName();

    } else {

        int reply2 = QMessageBox::question(
                    NULL, "Question", "Do you want an automatic segmentation?", QMessageBox::Yes, QMessageBox::No);

        if (reply2 == QMessageBox::Yes) {

            //Check for selection of image
            QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
            if (nodes.size() != 1) {
                QMessageBox::warning(NULL, "Attention", "Please select the CEMRA images from the Data Manager to segment!");
                return;
            }//_if

            //Test if this data item is an image
            QString mraPath;
            mitk::BaseData::Pointer data = nodes.at(0)->GetData();
            if (data) {
                mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
                if (image) {

                    this->BusyCursorOn();
                    mitk::ProgressBar::GetInstance()->AddStepsToDo(2);

                    //CNN prediction
                    mraPath = directory + mitk::IOUtil::GetDirectorySeparator() + "test.nii";
                    mitk::IOUtil::Save(image, mraPath.toStdString());
                    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
                    cmd->SetUseDockerContainers(true);
                    QString cnnPath = cmd->DockerCemrgNetPrediction(mraPath);
                    mitk::ProgressBar::GetInstance()->Progress();

                    //Clean prediction
                    using ImageTypeCHAR = itk::Image<short, 3>;
                    using ConnectedComponentImageFilterType = itk::ConnectedComponentImageFilter<ImageTypeCHAR, ImageTypeCHAR>;
                    using LabelShapeKeepNObjImgFilterType = itk::LabelShapeKeepNObjectsImageFilter<ImageTypeCHAR>;
                    ImageTypeCHAR::Pointer orgSegImage = ImageTypeCHAR::New();
                    mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(cnnPath.toStdString()), orgSegImage);
                    ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
                    connected->SetInput(orgSegImage);
                    connected->Update();
                    LabelShapeKeepNObjImgFilterType::Pointer lblShpKpNObjImgFltr = LabelShapeKeepNObjImgFilterType::New();
                    lblShpKpNObjImgFltr->SetInput(connected->GetOutput());
                    lblShpKpNObjImgFltr->SetBackgroundValue(0);
                    lblShpKpNObjImgFltr->SetNumberOfObjects(1);
                    lblShpKpNObjImgFltr->SetAttribute(LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
                    lblShpKpNObjImgFltr->Update();
                    mitk::Image::Pointer segImage = mitk::Image::New();
                    mitk::CastToMitkImage(lblShpKpNObjImgFltr->GetOutput(), segImage);
                    cnnPath = directory + mitk::IOUtil::GetDirectorySeparator() + "LA.nii";
                    mitk::IOUtil::Save(segImage, cnnPath.toStdString());
                    mitk::IOUtil::Load(cnnPath.toStdString(), *this->GetDataStorage());
                    remove(mraPath.toStdString().c_str());
                    fileName = "LA.nii";
                    QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
                    mitk::ProgressBar::GetInstance()->Progress();
                    this->BusyCursorOff();

                } else
                    QMessageBox::warning(NULL, "Attention", "Please select a CEMRA to segment!");
            }//_if_data

        } else {

            //Show the plugin
            this->GetSite()->GetPage()->ShowView("org.mitk.views.segmentation");

        }//_if_q2
    }//_if_q1
}

/**
 * @brief Relies on the interpolation mesh being added to DataManger
 * @param node Contains segmentation image to be saved
 */
void AtrialScarView::NodeAdded(const mitk::DataNode* node) {

    QString name = QString::fromStdString(node->GetName());
    if (name.endsWith("interpolation")) {

        mitk::DataStorage::SetOfObjects::ConstPointer nodes = this->GetDataStorage()->GetSources(node);
        mitk::DataStorage::SetOfObjects::ConstIterator itSeg = nodes->Begin();
        mitk::DataNode::Pointer segNode = itSeg->Value();

        if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

        //Find the selected node
        QString path;
        mitk::BaseData::Pointer data = segNode->GetData();
        if (data) {
            //Test if this data item is an image
            mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
            if (image) {

                bool ok;
                QString tmpFileName = fileName;
                fileName = QInputDialog::getText(
                            NULL, tr("Save Segmentation As"), tr("File Name:"), QLineEdit::Normal, fileName, &ok);
                if (ok && !fileName.isEmpty() && fileName.endsWith(".nii")) {
                    segNode->SetName(fileName.left(fileName.lastIndexOf(QChar('.'))).toStdString());
                    path = directory + mitk::IOUtil::GetDirectorySeparator() + fileName;
                    mitk::IOUtil::Save(image, path.toStdString());
                } else {
                    fileName = tmpFileName;
                    QMessageBox::warning(NULL, "Attention", "Please type a file name with the right extension (i.e. .nii)!");
                    return;
                }//_fileName

            } else {
                QMessageBox::warning(NULL, "Attention", "Segmentation image could not be saved automatically from the Data Manager!");
                return;
            }//_image
        } else
            return;
        this->GetDataStorage()->Remove(node);

    } else if (name.endsWith("PVeinsCroppedImage")) {

        //Listen to image crop change
        fileName = "PVeinsCroppedImage.nii";

    }//_if
}

void AtrialScarView::Registration() {

    //Toggle visibility of buttons
    if (m_Controls.button_x_1->isVisible()) {
        m_Controls.button_x_1->setVisible(false);
        m_Controls.button_x_2->setVisible(false);
    } else {
        m_Controls.button_x_1->setVisible(true);
        m_Controls.button_x_2->setVisible(true);
    }
}

void AtrialScarView::Register() {

    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(
                NULL, "Image Registration",
                "Have you completed steps 1 to 4 before using this feature?", QMessageBox::Yes|QMessageBox::No);
    if (reply == QMessageBox::No)
        return;

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 2) {
        MITK_INFO << ("Selection size:" + QString::number(nodes.size())).toStdString();
        QMessageBox::warning(NULL, "Attention", "Please select both LGE and CEMRA images from the Data Manager to register!");
        return;
    }

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    //Sort the two images
    QString lge, mra;
    mitk::DataNode::Pointer _1st = nodes.at(0);
    mitk::DataNode::Pointer _2nd = nodes.at(1);
    mitk::BaseData::Pointer dat1 = _1st->GetData();
    mitk::BaseData::Pointer dat2 = _2nd->GetData();

    if (dat1 && dat2) {
        //Test if this data items are image
        mitk::Image::Pointer image1 = dynamic_cast<mitk::Image*>(dat1.GetPointer());
        mitk::Image::Pointer image2 = dynamic_cast<mitk::Image*>(dat2.GetPointer());
        if (image1 && image2) {

            if (_1st->GetName().find("LGE") != _1st->GetName().npos && _2nd->GetName().find("MRA") != _2nd->GetName().npos) {

                lge = QString::fromStdString(_1st->GetName());
                mra = QString::fromStdString(_2nd->GetName());

            } else if (_1st->GetName().find("MRA") != _1st->GetName().npos && _2nd->GetName().find("LGE") != _2nd->GetName().npos) {

                lge = QString::fromStdString(_2nd->GetName());
                mra = QString::fromStdString(_1st->GetName());

            } else {

                QMessageBox::warning(NULL, "Attention", "Please select both LGE and CEMRA images from the Data Manager!");
                return;

            }//_if

            //Commandline call
            this->BusyCursorOn();
            mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
            std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
            cmd->SetUseDockerContainers(_useDockerInPlugin);
            cmd->ExecuteRegistration(directory, lge, mra);
            QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
            mitk::ProgressBar::GetInstance()->Progress();
            this->BusyCursorOff();

        } else {
            QMessageBox::warning(NULL, "Attention", "Please select both LGE and CEMRA images from the Data Manager!");
            return;
        }//_image
    } else
        return;
}

void AtrialScarView::Transform() {

    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(
                NULL, "Image Transformation",
                "Have you completed image registration before using this step?", QMessageBox::Yes|QMessageBox::No);
    if (reply == QMessageBox::No)
        return;

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1){
        MITK_WARN << ("[Transform] Problem with selection. Selection size: " + QString::number(nodes.size())).toStdString();
        QMessageBox::warning(NULL, "Attention", "Please select the corresponding segmentation to transform!");
        return;
    }

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    //Find the selected node
    QString path, pathTemp;
    mitk::DataNode::Pointer segNode = nodes.at(0);
    mitk::BaseData::Pointer data = segNode->GetData();
    if (data) {
        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {

            //Check seg node name
            QString segNodeTest = QString::fromStdString(segNode->GetName()) + ".nii";
            if (!fileName.contains(segNodeTest, Qt::CaseSensitive)) {
                MITK_WARN << ("[Transform] Problem with filename: " + fileName).toStdString();
                MITK_INFO << "[...][warning] segNode: " + segNode->GetName();
                QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
                return;
            }//_if

            bool ok;
            QString regFileName;
            regFileName = fileName.left(fileName.lastIndexOf(QChar('.'))) + "-reg.nii";
            regFileName = QInputDialog::getText(NULL, tr("Save Registration As"), tr("File Name:"), QLineEdit::Normal, regFileName, &ok);
            if (ok && !regFileName.isEmpty() && regFileName.endsWith(".nii")) {

                pathTemp = directory + mitk::IOUtil::GetDirectorySeparator() + "temp.nii";
                mitk::IOUtil::Save(image, pathTemp.toStdString());

                //Commandline call
                this->BusyCursorOn();
                mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
                std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
                cmd->SetUseDockerContainers(_useDockerInPlugin);
                cmd->ExecuteTransformation(directory, pathTemp.right(8), regFileName);
                QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
                mitk::ProgressBar::GetInstance()->Progress();
                this->BusyCursorOff();

                //Load the new segementation
                path = directory + mitk::IOUtil::GetDirectorySeparator() + regFileName;
                mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
                remove(pathTemp.toStdString().c_str());

                //Clear data manager
                mitk::DataStorage::SetOfObjects::ConstPointer sob = this->GetDataStorage()->GetAll();
                for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt)
                    if (nodeIt->Value()->GetName().find("MRA") != nodeIt->Value()->GetName().npos)
                        this->GetDataStorage()->Remove(nodeIt->Value());
                this->GetDataStorage()->Remove(segNode);
                fileName = regFileName;

            } else {
                QMessageBox::warning(NULL, "Attention", "Please type a file name with the right extension (i.e. .nii)!");
                return;
            }//_fileName

        } else {
            QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
            return;
        }//_image
    } else
        return;
}

void AtrialScarView::ClipPVeins() {

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    //Show the plugin
    this->GetSite()->GetPage()->ResetPerspective();
    AtrialScarClipperView::SetDirectoryFile(directory, fileName);
    this->GetSite()->GetPage()->ShowView("org.mitk.views.scarclipper");
}

void AtrialScarView::ExtraCalcs() {
    MITK_INFO << "[INFO] Calling view org.mitk.views.scarcalculations" << std::endl;
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    //Show the plugin
    this->GetSite()->GetPage()->ResetPerspective();
    ScarCalculationsView::SetCalculationsPaths(directory);
    if (!ScarCalculationsView::CheckForRequiredFiles()) {
        QMessageBox::warning(NULL, "Attention - Required files missing",
                             "Check all the files in the PRE/ANALYSIS/ and POST/ANALYSIS folders are correctly named.");
        directory = QString();
        return;
    }
    ScarCalculationsView::GetInputsFromFile();
    this->GetSite()->GetPage()->ShowView("org.mitk.views.scarcalculations");
}

void AtrialScarView::CreateSurf() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation to create a surface!");
        return;
    }

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    //Find the selected node
    QString path, pathTemp;
    mitk::DataNode::Pointer segNode = nodes.at(0);
    mitk::BaseData::Pointer data = segNode->GetData();
    if (data) {
        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {

            //Check seg node name
            if (segNode->GetName().compare(fileName.left(fileName.lastIndexOf(QChar('.'))).toStdString()) != 0) {
                QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
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

                bool ok1, ok2, ok3, ok4;
                int iter = m_UIMeshing.lineEdit_1->text().toInt(&ok1);
                float th = m_UIMeshing.lineEdit_2->text().toFloat(&ok2);
                int blur = m_UIMeshing.lineEdit_3->text().toInt(&ok3);
                int smth = m_UIMeshing.lineEdit_4->text().toInt(&ok4);

                //Set default values
                if (!ok1 || !ok2 || !ok3 || !ok4)
                    QMessageBox::warning(NULL, "Attention", "Reverting to default parameters!");
                if (!ok1) iter = 1;
                if (!ok2) th   = 0.5;
                if (!ok3) blur = 0;
                if (!ok4) smth = 10;
                //_if

                this->BusyCursorOn();
                pathTemp = directory + mitk::IOUtil::GetDirectorySeparator() + "temp.nii";
                mitk::IOUtil::Save(image, pathTemp.toStdString());
                mitk::ProgressBar::GetInstance()->AddStepsToDo(3);
                std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
                cmd->SetUseDockerContainers(_useDockerInPlugin);
                path = cmd->ExecuteSurf(directory, pathTemp, "close", iter, th, blur, smth);
                QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
                this->BusyCursorOff();

                //Add the mesh to storage
                std::string meshName = segNode->GetName() + "-Mesh";
                CemrgCommonUtils::AddToStorage(
                            CemrgCommonUtils::LoadVTKMesh(path.toStdString()), meshName, this->GetDataStorage());
                inputs->deleteLater();
                remove(pathTemp.toStdString().c_str());

            } else if (dialogCode == QDialog::Rejected) {
                inputs->close();
                inputs->deleteLater();
            }//_if

        } else {
            QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
            return;
        }//_image
    } else
        return;
}

void AtrialScarView::ClipperMV() {

    //Toggle visibility of buttons
    if (m_Controls.button_z_1->isVisible()) {
        m_Controls.button_z_1->setVisible(false);
        m_Controls.button_z_2->setVisible(false);
        return;
    } else {
        m_Controls.button_z_1->setVisible(true);
        m_Controls.button_z_2->setVisible(true);
    }//_if
}

void AtrialScarView::SelectLandmarks() {

    if (m_Controls.button_z_1->text() == QString::fromStdString("Select Landmarks")) {

        //Show the plugin
        QMessageBox::information(NULL, "Attention", "Please select 3 points around the mitral valve!");
        this->GetSite()->GetPage()->ShowView("org.mitk.views.pointsetinteraction");
        m_Controls.button_z_1->setText("Display Clipper");

    } else {

        //Check for selection of points
        QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
        if (nodes.empty()) {
            QMessageBox::warning(NULL, "Attention", "Please select the pointsets from the Data Manager to clip the mitral valve!");
            return;
        }//_if

        //Check selection type
        mitk::DataNode::Pointer landMarks = nodes.front();
        mitk::PointSet::Pointer pointSet = dynamic_cast<mitk::PointSet*>(landMarks->GetData());
        if (!pointSet || pointSet->GetSize() != 3) {
            QMessageBox::warning(NULL, "Attention", "Please select landmarks with 3 points from the Data Manager to continue!");
            return;
        }//_if

        //If the path was chosen incorrectly returns
        if (!RequestProjectDirectoryFromUser()) return;

        //Reset the button
        m_Controls.button_z_1->setText("Select Landmarks");

        //Retrieve mean and distance of 3 points
        double x_c = 0;
        double y_c = 0;
        double z_c = 0;
        for(int i=0; i<pointSet->GetSize(); i++) {
            x_c = x_c + pointSet->GetPoint(i).GetElement(0);
            y_c = y_c + pointSet->GetPoint(i).GetElement(1);
            z_c = z_c + pointSet->GetPoint(i).GetElement(2);
        }//_for
        x_c /= pointSet->GetSize();
        y_c /= pointSet->GetSize();
        z_c /= pointSet->GetSize();
        //double distance[pointSet->GetSize()];
        double * distance = new double[pointSet->GetSize()];
        for(int i=0; i<pointSet->GetSize(); i++) {
            double x_d = pointSet->GetPoint(i).GetElement(0) - x_c;
            double y_d = pointSet->GetPoint(i).GetElement(1) - y_c;
            double z_d = pointSet->GetPoint(i).GetElement(2) - z_c;
            distance[i] = sqrt(pow(x_d,2) + pow(y_d,2) + pow(z_d,2));
        }//_for
        double radius = *std::max_element(distance, distance + pointSet->GetSize());

        delete[] distance;
        //Create the clipper geometry
        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(x_c, y_c, z_c);
        sphereSource->SetRadius(radius);
        sphereSource->SetPhiResolution(40);
        sphereSource->SetThetaResolution(40);
        sphereSource->Update();
        mitk::Surface::Pointer mvClipper = mitk::Surface::New();
        mvClipper->SetVtkPolyData(sphereSource->GetOutput());

        //Adjust the data storage
        mitk::DataStorage::SetOfObjects::ConstPointer sob = this->GetDataStorage()->GetAll();
        for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt)
            if (nodeIt->Value()->GetName().find("MVClipper") != nodeIt->Value()->GetName().npos)
                this->GetDataStorage()->Remove(nodeIt->Value());
        CemrgCommonUtils::AddToStorage(mvClipper, "MVClipper", this->GetDataStorage(), false);
        sob = this->GetDataStorage()->GetAll();
        for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt) {
            if (nodeIt->Value()->GetName().find("MVClipper") != nodeIt->Value()->GetName().npos) {
                nodeIt->Value()->SetProperty("opacity", mitk::FloatProperty::New(0.4));
                nodeIt->Value()->SetProperty("color", mitk::ColorProperty::New(1.0, 0.0, 0.0));
            }//_if
        }//_for

    }//_if
}

void AtrialScarView::ClipMitralValve() {

    //Check for selection of points
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(NULL, "Attention", "Please select the pointsets from the Data Manager to clip the mitral valve!");
        this->GetSite()->GetPage()->ResetPerspective();
        return;
    }//_if

    //Check selection type
    mitk::DataNode::Pointer landMarks = nodes.front();
    mitk::PointSet::Pointer pointSet = dynamic_cast<mitk::PointSet*>(landMarks->GetData());
    if (!pointSet || pointSet->GetSize() != 3) {
        QMessageBox::warning(NULL, "Attention", "Please select landmarks with 3 points from the Data Manager to continue!");
        return;
    }//_if

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    //Read in and copy
    QString path = directory + mitk::IOUtil::GetDirectorySeparator() + "segmentation.vtk";
    mitk::Surface::Pointer surface = CemrgCommonUtils::LoadVTKMesh(path.toStdString());
    if (surface->GetVtkPolyData() == NULL) {
        QMessageBox::critical(NULL, "Attention", "No mesh was found in the project directory!");
        return;
    }//_if
    QString orgP = path.left(path.lastIndexOf(QChar('.'))) + "-Original.vtk";
    mitk::IOUtil::Save(mitk::IOUtil::Load<mitk::Surface>(path.toStdString()), orgP.toStdString());

    /*
     * Producibility Test
     **/
    QString prodPath = directory + mitk::IOUtil::GetDirectorySeparator();
    mitk::IOUtil::Save(pointSet, (prodPath + "prodMVCLandmarks.mps").toStdString());
    /*
     * End Test
     **/

    this->BusyCursorOn();
    std::unique_ptr<CemrgScar3D> scarObj = std::unique_ptr<CemrgScar3D>(new CemrgScar3D());
    surface = scarObj->ClipMesh3D(surface, pointSet);
    this->BusyCursorOff();

    //Check to remove the previous mesh node
    mitk::DataStorage::SetOfObjects::ConstPointer sob = this->GetDataStorage()->GetAll();
    for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt) {
        if (nodeIt->Value()->GetName().find("-Mesh") != nodeIt->Value()->GetName().npos)
            this->GetDataStorage()->Remove(nodeIt->Value());
        if (nodeIt->Value()->GetName().find("MVClipper") != nodeIt->Value()->GetName().npos)
            this->GetDataStorage()->Remove(nodeIt->Value());
    }//_for
    CemrgCommonUtils::AddToStorage(surface, "MVClipped-Mesh", this->GetDataStorage(), false);

    //Reverse coordination of surface for writing MIRTK style
    mitk::Surface::Pointer surfCloned = surface->Clone();
    vtkSmartPointer<vtkPolyData> pd = surfCloned->GetVtkPolyData();
    for (int i=0; i<pd->GetNumberOfPoints(); i++) {
        double* point = pd->GetPoint(i);
        point[0] = -point[0];
        point[1] = -point[1];
        pd->GetPoints()->SetPoint(i, point);
    }//_for
    mitk::IOUtil::Save(surfCloned, path.toStdString());
}

void AtrialScarView::ScarMap() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(NULL, "Attention", "Please select the LGE image from the Data Manager to calculate the scar map!");
        return;
    }
    std::size_t found = nodes.at(0)->GetName().find("LGE");
    bool check = found != std::string::npos ? false : true;
    if (check) {
        QMessageBox::warning(NULL, "Attention", "Please select the LGE image from the Data Manager to calculate the scar map!");
        return;
    }

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    //Check for mesh in the project directory
    try {
        QString path = directory + mitk::IOUtil::GetDirectorySeparator() + "segmentation.vtk";
        mitk::IOUtil::Load<mitk::Surface>(path.toStdString());
    } catch (...) {
        QMessageBox::critical(NULL, "Attention", "No mesh was found in the project directory!");
        return;
    }//_try

    //Find the selected node
    mitk::DataNode::Pointer imgNode = nodes.at(0);
    mitk::BaseData::Pointer data = imgNode->GetData();
    if (data) {

        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image.IsNotNull()) {

            scar = std::unique_ptr<CemrgScar3D>(new CemrgScar3D());
            if (scar) {

                //Ask for user input to set the parameters
                QDialog* inputs = new QDialog(0,0);
                m_UIScar.setupUi(inputs);
                connect(m_UIScar.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
                connect(m_UIScar.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
                int dialogCode = inputs->exec();

                //Act on dialog return code
                if (dialogCode == QDialog::Accepted) {

                    bool ok1, ok2;
                    int minStep = m_UIScar.lineEdit_1->text().toInt(&ok1);
                    int maxStep = m_UIScar.lineEdit_2->text().toInt(&ok2);
                    int methodType = m_UIScar.radioButton_1->isChecked() ? 2 : 1;
                    QString meType = m_UIScar.radioButton_1->isChecked() ? "Max" : "Mean";

                    //Set default values
                    if (!ok1 || !ok2)
                        QMessageBox::warning(NULL, "Attention", "Reverting to default parameters!");
                    if (!ok1) minStep =-1;
                    if (!ok2) maxStep = 3;
                    //_if

                    /*
                     * Producibility Test
                     **/
                    QString prodPath = directory + mitk::IOUtil::GetDirectorySeparator();
                    ofstream prodFile1;
                    prodFile1.open((prodPath + "prodScarMapInputs.txt").toStdString());
                    prodFile1 << minStep << "\n";
                    prodFile1 << maxStep << "\n";
                    prodFile1 << meType.toStdString() << "\n";
                    prodFile1.close();
                    /*
                     * End Test
                     **/

                    this->BusyCursorOn();
                    mitk::ProgressBar::GetInstance()->AddStepsToDo(1);

                    scar->SetMinStep(minStep);
                    scar->SetMaxStep(maxStep);
                    scar->SetMethodType(methodType);
                    mitk::Image::Pointer scarSegImg;
                    try {
                        QString path = directory + mitk::IOUtil::GetDirectorySeparator() + fileName;
                        scarSegImg = mitk::IOUtil::Load<mitk::Image>(path.toStdString());
                    } catch (...) {
                        QMessageBox::critical(NULL, "Attention", "The loaded or created segmentation was not found!");
                        mitk::ProgressBar::GetInstance()->Progress();
                        this->BusyCursorOff();
                        return;
                    }//_try

                    //Resample to fit LGE
                    typedef itk::Image<short, 3> ImageType;
                    itk::ResampleImageFilter<ImageType, ImageType>::Pointer resampleFilter;
                    ImageType::Pointer scarSegITK = ImageType::New();
                    CastToItkImage(scarSegImg, scarSegITK);
                    ImageType::Pointer lgeITK = ImageType::New();
                    CastToItkImage(image, lgeITK);
                    resampleFilter = itk::ResampleImageFilter<ImageType, ImageType >::New();
                    resampleFilter->SetInput(scarSegITK);
                    resampleFilter->SetReferenceImage(lgeITK);
                    resampleFilter->SetUseReferenceImage(true);
                    resampleFilter->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<ImageType>::New());
                    resampleFilter->SetDefaultPixelValue(0);
                    resampleFilter->UpdateLargestPossibleRegion();
                    scarSegImg = mitk::ImportItkImage(resampleFilter->GetOutput());

                    //Update datamanger
                    std::string nodeSegImgName = fileName.left(fileName.lastIndexOf(QChar('.'))).toStdString();
                    mitk::DataStorage::SetOfObjects::ConstPointer sob = this->GetDataStorage()->GetAll();
                    for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt)
                        if (nodeIt->Value()->GetName().find(nodeSegImgName) != nodeIt->Value()->GetName().npos)
                            this->GetDataStorage()->Remove(nodeIt->Value());
                    QString savePath = directory + mitk::IOUtil::GetDirectorySeparator() + fileName;
                    mitk::IOUtil::Save(scarSegImg, savePath.toStdString());

                    //Projection
                    scar->SetScarSegImage(scarSegImg);
                    mitk::Surface::Pointer shell = scar->Scar3D(directory.toStdString(), image);
                    mitk::DataNode::Pointer node = CemrgCommonUtils::AddToStorage(shell, (meType + "Scar3D").toStdString(), this->GetDataStorage());

                    MITK_INFO << "Saving debug scar map labels.";
                    scar->SaveScarDebugImage(meType + "_debugSCAR.nii", directory);

                    //Check to remove the previous mesh node
                    sob = this->GetDataStorage()->GetAll();
                    for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt)
                        if (nodeIt->Value()->GetName().find("-Mesh") != nodeIt->Value()->GetName().npos)
                            this->GetDataStorage()->Remove(nodeIt->Value());

                    //LUT setup
                    mitk::LookupTable::Pointer scarLut = mitk::LookupTable::New();
                    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
                    lut->SetTableRange(0, scar->GetMaxScalar());
                    lut->SetHueRange(0.33, 0.0);
                    lut->Build();
                    lut->SetTableValue(0, 0.0, 0.0, 0.0, 1.0);
                    scarLut->SetVtkLookupTable(lut);
                    node->SetProperty("LookupTable", mitk::LookupTableProperty::New(scarLut));
                    node->SetProperty("scalar visibility", mitk::BoolProperty::New(true));
                    node->SetProperty("scalar mode", mitk::VtkScalarModeProperty::New(2));
                    node->SetFloatProperty("ScalarsRangeMinimum", 0);
                    node->SetFloatProperty("ScalarsRangeMaximum", scar->GetMaxScalar());

                    //Save the vtk mesh
                    QString name(imgNode->GetName().c_str());
                    name = name.right(name.length() - name.lastIndexOf("-") - 1);
                    QString path = directory + mitk::IOUtil::GetDirectorySeparator() + name + "-" + meType + "Scar.vtk";
                    mitk::IOUtil::Save(shell, path.toStdString());

                    mitk::ProgressBar::GetInstance()->Progress();
                    this->BusyCursorOff();
                    inputs->deleteLater();

                } else if (dialogCode == QDialog::Rejected) {
                    inputs->close();
                    inputs->deleteLater();
                }//_if
            }//_scar
        } else
            return;
    } else
        return;
}

void AtrialScarView::ScarDebug() {

    MITK_INFO << "Loading: " + debugSCARname.toStdString();

    if (debugSCARname.isEmpty()) {
        MITK_INFO << "File: " + debugSCARname.toStdString() + " NOT FOUND";
        return;
    }
    mitk::IOUtil::Load(debugSCARname.toStdString(), *this->GetDataStorage());
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
    m_Controls.button_deb->setVisible(false);
    //Restore image name
    //char sep = mitk::IOUtil::GetDirectorySeparator();
    //fileName = path.mid(path.lastIndexOf(sep) + 1);
}

void AtrialScarView::Threshold() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(NULL, "Attention", "Please select the LGE image from the Data Manager to quantify the scar!");
        return;
    }
    std::size_t found = nodes.at(0)->GetName().find("LGE");
    bool check = found != std::string::npos ? false : true;
    if (check) {
        QMessageBox::warning(NULL, "Attention", "Please select the LGE image from the Data Manager to quantify the scar!");
        return;
    }

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    //Find the selected node
    double mean = 0.0, stdv = 0.0;
    mitk::DataNode::Pointer imgNode = nodes.at(0);
    mitk::BaseData::Pointer data = imgNode->GetData();
    if (data) {

        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image.IsNotNull()) {

            //Convert images to right type
            mitk::Image::Pointer roi;
            mitk::Image::Pointer lgeImage = mitk::Image::New();
            itk::Image<float,3>::Pointer itkImage = itk::Image<float,3>::New();
            mitk::CastToItkImage(image, itkImage);
            mitk::CastToMitkImage(itkImage, lgeImage);
            try {
                QString path = directory + mitk::IOUtil::GetDirectorySeparator() + fileName;
                roi = mitk::IOUtil::Load<mitk::Image>(path.toStdString());
            } catch (...) {
                QMessageBox::critical(NULL, "Attention", "The loaded or created segmentation was not found!");
                return;
            }//_try
            mitk::Image::Pointer roiImage = mitk::Image::New();
            itk::Image<float,3>::Pointer roiItkImage = itk::Image<float,3>::New();
            mitk::CastToItkImage(roi, roiItkImage);
            if (scar) {

                //Erosion of bloodpool
                typedef itk::Image<float, 3> ImageType;
                typedef itk::BinaryCrossStructuringElement<ImageType::PixelType, 3> CrossType;
                typedef itk::GrayscaleErodeImageFilter<ImageType, ImageType, CrossType> ErosionFilterType;
                bool ok;
                int vxls = QInputDialog::getInt(NULL, tr("Bloodpool Erosion"), tr("Voxels:"), 3, 1, 20, 1, &ok);
                if (!ok) {
                    QMessageBox::warning(NULL, "Attention", "Enter a correct value!");
                    return;
                }//_if
                CrossType binaryCross;
                binaryCross.SetRadius(vxls);
                binaryCross.CreateStructuringElement();
                ErosionFilterType::Pointer erosionFilter = ErosionFilterType::New();
                erosionFilter->SetInput(roiItkImage);
                erosionFilter->SetKernel(binaryCross);
                erosionFilter->UpdateLargestPossibleRegion();
                roiImage = mitk::ImportItkImage(erosionFilter->GetOutput())->Clone();
                CemrgCommonUtils::AddToStorage(roiImage, "Eroded ROI", this->GetDataStorage());

                //Calculate mean, std of ROI
                bool success = scar->CalculateMeanStd(lgeImage, roiImage, mean, stdv);
                if (!success)
                    return;

            } else {
                QMessageBox::warning(NULL, "Attention", "The scar map from the previous step has not been generated!");
                return;
            }//_if_scar
        }//_image
    }//_data

    //Ask for user input to set the parameters
    QDialog* inputs = new QDialog(0,0);
    m_UISQuant.setupUi(inputs);
    connect(m_UISQuant.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UISQuant.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    int dialogCode = inputs->exec();

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {

        bool ok;
        double value = m_UISQuant.lineEdit->text().toDouble(&ok);
        int methodType = m_UISQuant.radioButton_1->isChecked() ? 1 : 2;

        //Set default values
        if (!ok) {
            QMessageBox::warning(NULL, "Attention", "Please enter a valid value!");
            return;
        }//_ok

        this->BusyCursorOn();
        mitk::ProgressBar::GetInstance()->AddStepsToDo(1);

        double percentage = 0;
        double thresh = (methodType == 1) ? mean*value : mean+value*stdv;

        /*
         * Producibility Test
         **/
        QString prodPath = directory + mitk::IOUtil::GetDirectorySeparator();
        ofstream prodFile1;
        prodFile1.open((prodPath + "prodThresholds.txt").toStdString());
        prodFile1 << value << "\n";
        prodFile1 << methodType << "\n";
        prodFile1 << mean << "\n";
        prodFile1 << stdv << "\n";
        prodFile1 << thresh << "\n";
        prodFile1.close();
        /*
         * End Test
         **/

        if (scar) percentage = scar->Thresholding(thresh);
        std::ostringstream os;
        os << std::fixed << std::setprecision(2) << percentage;
        QString message = "The percentage scar is " + QString::fromStdString(os.str()) + "% of total segmented volume.";
        QMessageBox::information(NULL, "Scar Quantification", message);

        mitk::ProgressBar::GetInstance()->Progress();
        this->BusyCursorOff();
        inputs->deleteLater();

    } else if (dialogCode == QDialog::Rejected) {
        inputs->close();
        inputs->deleteLater();
    }//_if
}

void AtrialScarView::Sphericity() {

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    //Read in the mesh
    QString path = directory + mitk::IOUtil::GetDirectorySeparator() + "segmentation.vtk";
    mitk::Surface::Pointer surface = CemrgCommonUtils::LoadVTKMesh(path.toStdString());

    double result = 0;
    std::unique_ptr<CemrgMeasure> Sphericity = std::unique_ptr<CemrgMeasure>(new CemrgMeasure());
    if (surface->GetVtkPolyData() != NULL)
        result = Sphericity->GetSphericity(surface->GetVtkPolyData());
    else {
        QMessageBox::warning(NULL, "Attention", "No mesh was found in the project directory!");
        return;
    }//_if

    QString message = "The sphericity is " + QString::number(result) + "%.";
    QMessageBox::information(NULL, "Sphericity", message);
}

void AtrialScarView::ResetMain() {

    Reset(true);
    m_Controls.button_3->setVisible(true);
    m_Controls.button_4->setVisible(false);
    m_Controls.button_5->setVisible(false);
    m_Controls.button_6->setVisible(false);
    m_Controls.button_7->setVisible(false);
    m_Controls.button_x->setVisible(false);
    m_Controls.button_y->setVisible(false);
    m_Controls.button_z->setVisible(false);
    m_Controls.button_s->setVisible(false);
    m_Controls.button_2_1->setVisible(false);
    m_Controls.button_x_1->setVisible(false);
    m_Controls.button_x_2->setVisible(false);
    m_Controls.button_z_1->setVisible(false);
    m_Controls.button_z_2->setVisible(false);
    m_Controls.button_deb->setVisible(false);
}

void AtrialScarView::Reset(bool allItems) {

    try {

        ctkPluginContext* context = mitk::kcl_cemrgapp_scar_Activator::getContext();
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

        //Remove everything or keep images
        mitk::DataStorage::SetOfObjects::ConstPointer nodesToRemove = dataStorage->GetAll();
        if (allItems) {
            dataStorage->Remove(nodesToRemove);
        } else {
            mitk::DataStorage::SetOfObjects::ConstIterator iter = nodesToRemove->Begin();
            mitk::DataStorage::SetOfObjects::ConstIterator iterEnd = nodesToRemove->End();
            for (; iter != iterEnd; ++iter) {
                QString name = QString::fromStdString(iter->Value()->GetName());
                if (name.startsWith("dcm-"))
                    continue;
                dataStorage->Remove(iter->Value());
            }//_for
        }//_if

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
    fileName = ".nii";
    directory.clear();
    m_Controls.button_z_1->setText("Select Landmarks");
    this->GetSite()->GetPage()->ResetPerspective();
}

bool AtrialScarView::RequestProjectDirectoryFromUser() {

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
    }//_if

    return succesfulAssignment;
}
